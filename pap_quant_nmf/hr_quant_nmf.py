import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from scipy import stats
import seaborn as sns
from sklearn.metrics import r2_score
import matplotlib.ticker as mticker
from matplotlib.offsetbox import AnchoredText
import os
import warnings
import pickle
warnings.filterwarnings("ignore")


# an auxiliary function to compute the projection of a point on a line
# p1 and p2 define a line and the function aims to find a point p4 on the line that is closest to p3
def get_closet_point_in_line(p1, p2, p3):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    dx, dy = x2 - x1, y2 - y1
    det = dx * dx + dy * dy
    a = (dy * (y3 - y1) + dx * (x3 - x1)) / det
    return x1 + a * dx, y1 + a * dy


# parse command line argument
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename',
                    action='store',
                    dest='filename',
                    type=str,
                    help="pap data file name (csv or tsv files are preferred)")
parser.add_argument('-m', '--pretrained_nmf_model',
                    action='store',
                    dest='pretrained_nmf_model',
                    default=None,
                    type=str,
                    help="file name of pretrained nmf model (a pickle object)")
parser.add_argument('-mina', '--min_alpha',
                    action='store',
                    dest='min_alpha',
                    default=1e-4,
                    type=float,
                    help="minium value that multiplies the regularization terms")
parser.add_argument('-maxa', '--max_alpha',
                    action='store',
                    dest='max_alpha',
                    default=1e0,
                    type=float,
                    help="maximum value that multiplies the regularization terms")
parser.add_argument('-na', '--n_alphas',
                    action='store',
                    dest='n_alphas',
                    default=25,
                    type=int,
                    help="number of multipliers between min_alpha and max_alpha")
parser.add_argument('-tol', '--tolerance_deviate_from_optimum',
                    action='store',
                    dest='tol',
                    default=0.05,
                    type=float,
                    help="percent tolerance when alpha increases the reconstruction error)")
parser.add_argument('-o', '--output_directory',
                    action='store',
                    dest='out_dir',
                    help='output_directory (no trailing /)',
                    type=str,
                    default='output')
args = parser.parse_args()

# read PAP data in a text file
# input file should be a 2-dimensional table where rows are drug concentrations and columns are strains
if not args.filename:
    print('PAP file must be provided. Exit.\n')
    parser.print_help()
    exit()
if args.filename.endswith('.csv'):
    df_pap = pd.read_csv(args.filename, index_col=0)
elif args.filename.endswith('.tsv'):
    df_pap = pd.read_csv(args.filename, index_col=0, sep="\t")
else:
    warnings.warn("Unknown file extension. Use white space(s) as deliminator.")
    df_pap = pd.read_csv(args.filename, index_col=0, sep="\s+")

# make sure that the row names of df_pap are numerics (drug concentrations)
df_pap = df_pap.rename({idx: float(idx) for idx in df_pap.index})
df_pap.index.name = 'drug_conc'
df_pap = df_pap.reset_index().sort_values('drug_conc', ascending=True).set_index('drug_conc')

# load pretrained nmf model
if args.pretrained_nmf_model is not None:
    with open(args.pretrained_nmf_model, 'rb') as pickle_file:
        model = pickle.load(pickle_file)

# check min alpha
min_alpha = float(args.min_alpha)
if min_alpha <= 0.0:
    raise RuntimeError('min_alpha (%2.2f) must be positive.' % min_alpha)

# check max alpha
max_alpha = float(args.max_alpha)
if max_alpha < min_alpha:
    raise RuntimeError('max_alpha (%2.2f) must be no smaller than min alpha (%2.2f).' % (max_alpha, min_alpha))

# check n alphas
n_alphas = int(args.n_alphas)
if n_alphas <= 0:
    raise RuntimeError('n_alphas (%d) must be positive.' % n_alphas)

# check tolerance P2I
tol = float(args.tol)
if tol < 0.0:
    raise RuntimeError('tol (%2.2f) must be positive.' % tol)

# check if output folder exists. if not, create a new one
out_dir = args.out_dir
if not os.path.isdir(out_dir):
    os.system('mkdir %s' % out_dir)

# get detection limit
df_pap_stacked = df_pap.stack().reset_index()
df_pap_stacked.columns = ['drug_conc', 'strain','survival']
min_survival = df_pap_stacked[df_pap_stacked['survival'] > 0.0]['survival'].min()
dol = 10**np.floor(np.log10(min_survival))
shifted_max_log10_pap = 2 - np.log10(dol)
print('Min survival percentage = %2.2e. Detection of limit = %2.2e.' % (min_survival, dol))

# convert to log scale and 0 will be replaced by the detection of limit
with np.errstate(divide='ignore'):
    df_pap = np.log10(df_pap).replace(-np.inf, np.log10(dol)) - np.log10(dol)

# run non-negative matrix factorization
df_pap_t = df_pap.T
if args.pretrained_nmf_model is None:
    print("Running original NMF ...")
    success = False
    alphas = [0.0] + [10**x for x in np.linspace(np.log10(min_alpha), np.log10(max_alpha), n_alphas)]
    min_reconstruction_err = -1.0
    best_alpha = 0.0
    for alpha in alphas:
        model = NMF(n_components=2, init='nndsvdar', random_state=42, max_iter=100000, tol=1e-6, alpha_H=alpha, alpha_W=alpha)
        df_pap_W = pd.DataFrame(model.fit_transform(df_pap_t.values), columns=['Factor1', 'Factor2'], index=df_pap_t.index)
        df_pap_W.index.name = 'strain'
        df_pap_H = pd.DataFrame(model.components_, columns=df_pap.index, index=['Basis1', 'Basis2'])
        if model.n_iter_ == 100000:
            print('alpha = %2.2e, solution does not converge.' % alpha)
        else:
            print('alpha = %2.2e, reconstruction_error = %2.2f.' % (alpha, model.reconstruction_err_))
            if alpha == 0.0:
                min_reconstruction_err = model.reconstruction_err_
            else:
                if model.reconstruction_err_ > min_reconstruction_err * (1+tol):
                    # retrain the model using best alpha
                    model = NMF(n_components=2, init='nndsvdar', random_state=42, max_iter=100000, tol=1e-6, alpha_H=best_alpha, alpha_W=best_alpha)
                    df_pap_W = pd.DataFrame(model.fit_transform(df_pap_t.values), columns=['Factor1', 'Factor2'], index=df_pap_t.index)
                    df_pap_W.index.name = 'strain'
                    df_pap_H = pd.DataFrame(model.components_, columns=df_pap.index, index=['Basis1', 'Basis2'])

                    # save model as an object
                    success = True
                    with open('%s/best_model.pickle' % out_dir, 'wb') as handle:
                        pickle.dump(model, handle)
                    break
                else:
                    best_alpha = alpha
    if success == False:
        print('NMF fails. Exit.')
        exit()
else:
    print("Running NMF with precomputed basis ...")
    df_pap_W = pd.DataFrame(model.transform(df_pap_t.values), columns=['Factor1', 'Factor2'], index=df_pap_t.index)
    df_pap_W.index.name = 'strain'
    df_pap_H = pd.DataFrame(model.components_, columns=df_pap.index, index=['Basis1', 'Basis2'])

# normalize W and H matrices
multiplier = [shifted_max_log10_pap / v for v in list(df_pap_H[0])]
df_pap_H.iloc[0, :] = df_pap_H.iloc[0, :] * multiplier[0]
df_pap_H.iloc[1, :] = df_pap_H.iloc[1, :] * multiplier[1]
df_pap_W.iloc[:, 0] = df_pap_W.iloc[:, 0] / multiplier[0]
df_pap_W.iloc[:, 1] = df_pap_W.iloc[:, 1] / multiplier[1]

# determine which feature represents the S/HR state
basis1_feature = 'S'
basis2_feature = 'HR'
if np.mean(df_pap_H.loc['Basis1', :]) > np.mean(df_pap_H.loc['Basis2', :]):
    basis1_feature = 'HR'
    basis2_feature = 'S'

# compute hetero-resistance index by projecting all weights (W) on x+y=1
res = []
for strain in df_pap_W.index:
    projection = get_closet_point_in_line((0, 1), (1, 0),
                                          (df_pap_W.loc[strain, 'Factor1'], df_pap_W.loc[strain, 'Factor2']))
    if basis1_feature == 'HR':
        hr_index = projection[0]
    else:
        hr_index = projection[1]
    if hr_index < 0:
        warnings.warn("Hetero-resistance index < 0 for %s. its value (%2.2f) is round to 0.0." % (strain, hr_index))
        hr_index = 0
    if hr_index > 1:
        warnings.warn("Hetero-resistance index > 1 for %s. its value (%2.2f) is round to 1.0" % (strain, hr_index))
        hr_index = 1
    res.append(hr_index)
df_hr = pd.DataFrame(res, columns=['HR_INDEX'], index=df_pap_W.index).sort_values('HR_INDEX')

# plot
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(18, 12))
ax00 = plt.subplot(2, 3, 1)
ax01 = plt.subplot(2, 3, 2)
ax02 = plt.subplot(2, 3, 3)
ax10 = plt.subplot(2, 3, 4)
ax1_12 = plt.subplot2grid((2, 3), (1, 1), colspan=2)

# plot PAP curves
for strain in df_pap.columns:
    _ = ax00.plot(df_pap.index, df_pap[strain], linewidth=1, color='Gray')
_ = ax00.set_xlabel('Drug concentration', fontsize=12)
_ = ax00.set_ylabel('$log_{10}(Survival~\%)$', fontsize=12)
_ = ax00.set_xticks(list(df_pap.index))
_ = ax00.set_ylim([-0.5, shifted_max_log10_pap + 0.5])
_ = ax00.set_yticks(np.arange(0, shifted_max_log10_pap + 1))
_ = ax00.set_yticklabels([y + np.log10(dol) for y in np.arange(0, shifted_max_log10_pap + 1)])

# plot features (H matrix)
_ = ax01.plot(df_pap.index, df_pap_H.loc['Basis1', :], 'ro-', linewidth=2, markersize=8,
              label='%s feature' % basis1_feature)
_ = ax01.plot(df_pap.index, df_pap_H.loc['Basis2', :], 'bo-', linewidth=2, markersize=8,
              label='%s feature' % basis2_feature)
_ = ax01.legend()
_ = ax01.set_xlabel('Drug concentration', fontsize=12)
_ = ax01.set_ylabel('$log_{10}(Survival ~\%)$', fontsize=12)
max_yint = np.ceil(np.max(df_pap_H.max().max()))
_ = ax01.set_ylim([-0.5, max_yint + 0.5])
_ = ax01.set_yticks(np.arange(0, max_yint + 1))
_ = ax01.set_yticklabels([y + np.log10(dol) for y in np.arange(0, max_yint + 1)])
_ = ax01.set_xticks(list(df_pap.index))

# plot correlation of weights (W matrix)
regression = stats.linregress(df_pap_W.Factor1, df_pap_W.Factor2)
_ = sns.scatterplot(x='Factor1', y='Factor2', data=df_pap_W, color='gray', marker='o', s=100, ax=ax02,
                    label='R2 = %2.2f' % (regression.rvalue ** 2))
_ = ax02.legend()
xdense = np.linspace(-0.1, 1.1, 100)
ydense = regression.intercept + regression.slope * xdense
_ = ax02.plot(xdense, ydense, 'k--')
_ = ax02.set_xlabel("%s feature weight" % basis1_feature, fontsize=12)
_ = ax02.set_ylabel('%s feature weight' % basis2_feature, fontsize=12)
_ = ax02.set_xlim([-0.1, 1.1])
_ = ax02.set_ylim([-0.1, 1.1])
_ = ax02.legend()

# plot reconstruction error
res = []
for strain, hr in zip(df_hr.index, df_hr.HR_INDEX):
    pap_obs = df_pap[strain]
    pap_pred = df_pap_W.loc[strain, 'Factor1'] * df_pap_H.loc['Basis1', :] + \
               df_pap_W.loc[strain, 'Factor2'] * df_pap_H.loc['Basis2', :]
    for conc, obs, pred in zip(list(df_pap.index), pap_obs, pap_pred):
        res.append([strain, conc, obs, pred])
df_recon = pd.DataFrame(res, columns=['strain', 'drug_conc', 'obs', 'pred'])
df_recon['drug_conc'] = df_recon['drug_conc'].astype(str)
r2 = r2_score(y_true=df_recon['obs'], y_pred=df_recon['pred'])
_ = sns.scatterplot(x='obs', y='pred', hue='drug_conc', palette='Blues', ax=ax10, data=df_recon, s=36)
_ = ax10.set_xlim([-0.5, shifted_max_log10_pap + 0.5])
_ = ax10.set_ylim([-0.5, shifted_max_log10_pap + 0.5])
_ = ax10.plot([-0.5, shifted_max_log10_pap + 0.5], [-0.5, shifted_max_log10_pap + 0.5], 'k--')
_ = ax10.set_xlabel('Observed $log_{10}(Survival ~\%)$', fontsize=12)
_ = ax10.set_xticks(np.arange(0, shifted_max_log10_pap + 1))
_ = ax10.set_xticklabels([x + np.log10(dol) for x in np.arange(0, shifted_max_log10_pap + 1)])
_ = ax10.set_ylabel('Reconstructed $log_{10}(Survival ~\%)$', fontsize=12)
_ = ax10.set_yticks(np.arange(0, shifted_max_log10_pap + 1))
_ = ax10.set_yticklabels([y + np.log10(dol) for y in np.arange(0, shifted_max_log10_pap + 1)])
_ = ax10.set_title("R-squared: %2.2f" % r2)
text_box = AnchoredText('R2: %2.2f' % r2, frameon=True, loc=4, pad=0.5)
plt.setp(text_box.patch, facecolor='white', alpha=0.5)
_ = ax10.add_artist(text_box)

# plot bars of hetero-resistance index (ranked from low to high)
_ = sns.barplot(x=df_hr.index,
                y=df_hr.HR_INDEX,
                order=df_hr.index,
                ax=ax1_12,
                palette='Reds')
for item in ax1_12.get_xticklabels():
    item.set_rotation(90)
_ = ax1_12.set_xlabel('strain', fontsize=12)
_ = ax1_12.set_ylabel('HR index', fontsize=12)
_ = ax1_12.set_ylim([0, 1])
for patch in ax1_12.patches:
    current_width = patch.get_width()
    patch.set_width(1)
    patch.set_x(patch.get_x() + current_width - 1)
    patch.set_edgecolor('k')
max_num_to_display = 60
if len(df_hr) > max_num_to_display:
    myLocator = mticker.MultipleLocator(np.ceil(len(df_hr) / max_num_to_display))
    _ = ax1_12.xaxis.set_major_locator(myLocator)

# save HR INDEX to file
df_hr.to_csv('%s/hr_quant.csv' % out_dir)

# save to svg and png files
fig.savefig("%s/hr_quant.png" % out_dir, format="png", dpi=600)

# uncomment the following lines if a svg format is needed
# plt.rcParams['svg.fonttype'] = 'none'
# fig.savefig("hr_quant.svg", format="svg", dpi=600)
