import pandas as pd
import math
import os
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from numpy.random import lognormal, exponential, uniform
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

#########################
# PAP decomposition model
#########################
def pap_simulation(drug_concs, *p):

    global time_length
    global detection_limit
    global use_one_population_model
    params = [x for x in p]

    pred_cfu_ratio = []
    if use_one_population_model:
        cfu_ratio = np.exp(-(params[0]*drug_concs**params[1]/(drug_concs**params[1]+params[2]**params[1]))*time_length)
    else:
        # swap parameter 1,2 and 4,5 to guarantee that params[0] describes frequency of susceptible population
        if params[1]<params[2]:
            tmp = params[1]
            params[1] = params[2]
            params[2] = tmp
        if params[4]>params[5]:
            tmp = params[4]
            params[4] = params[5]
            params[5] = tmp
        cfu_ratio = params[0]*np.exp(-(params[1]*drug_concs**params[3]/(drug_concs**params[3]+params[4]**params[3]))*time_length) \
                  + (1-params[0])*np.exp(-(params[2]*drug_concs**params[3]/(drug_concs**params[3]+params[5]**params[3]))*time_length)

    # replace actual cfu ratio with detection limit when it is below detection limit
    for r in cfu_ratio:
        if r < detection_limit:
            pred_cfu_ratio.append(np.log10(detection_limit))
        else:
            pred_cfu_ratio.append(np.log10(r))

    return pred_cfu_ratio

#########################
# READ configuration file
#########################
if os.path.exists("config.csv"):
    df_config = pd.read_csv("config.csv", index_col=0, header=None)
    if 'pap_file' in list(df_config.index):
        pap_file = list(df_config.loc['pap_file'])[0]
        df_pap = pd.read_csv(pap_file, index_col=0)/100 # first column is drug concentration
        drug_concs = list(df_pap.index)
    else:
        raise RuntimeError ("pap_file (file path to pap data) not found.")

    if 'detection_limit' in list(df_config.index):
        detection_limit = float(list(df_config.loc['detection_limit'])[0])
    else:
        raise RuntimeError ("detection_limit (minimum cfu ratio) not provided.")

    if 'time_length' in list(df_config.index):
        time_length = float(list(df_config.loc['time_length'])[0])
    else:
        raise RuntimeError ("time_length (time duration of PAP) not provided.")

    if 'mc_runs' in list(df_config.index):
        mc_runs = int(list(df_config.loc['mc_runs'])[0])
    else:
        raise RuntimeError ("mc_runs (number of random sampling of initial conditions) not provided.")

    if 'mic_to_plot' in list(df_config.index):
        mic_to_plot = [float(x) for x in list(df_config.loc['mic_to_plot'])[0].split(';')]
        if len(mic_to_plot) >20:
            raise RuntimeError ("too many mic (only 20 allowed).")
    else:
        raise RuntimeError ("mic_to_plot (potential MIC values of susceptible population) not provided.")

    # optional
    if 'to_include' in list(df_config.index):
        strains_to_include = list(df_config.loc['to_include'])[0].split(';')
        df_pap = df_pap[[col for col in df_pap.columns if col in strains_to_include]]

    if 'maxg_ub' in list(df_config.index):
        maxg_ub = float(list(df_config.loc['maxg_ub'])[0])
    else:
        maxg_ub = np.infty
else:
    raise RuntimeError ("config.csv not found.")

####################
# DECOMPOSE PAP data
####################

results = []
for strain in df_pap.columns:
    print("-----------------")
    print("-----------------")
    print("strain:", strain)

    #--------------------
    # remove missing data
    #--------------------
    obs_cfu_ratio = list(np.log10(df_pap[strain]))
    xdata = []
    ydata = []
    for x,y in zip(drug_concs, obs_cfu_ratio):
        if not math.isnan(y):
            xdata.append(x)
            ydata.append(y)

    #----------------------------
    # run one population model
    #----------------------------
    print("\nrunning one population model...")
    use_one_population_model=True
    deltag = lognormal(0, 1, size=mc_runs)
    hillc = uniform(1, 4, size=mc_runs)
    K50 = exponential(10, size=mc_runs)
    lb = [0.0, 1.0, 0.0]
    ub = [np.inf, np.inf, np.inf]
    mse = []
    popt_list = []
    _iter = 1
    for p1, p2, p3 in zip(deltag, hillc, K50):
        curr_params = [p1, p2, p3]
        try:
            popt, pcov, infodict, mesg, ier = curve_fit(
                f=pap_simulation,
                xdata=xdata,
                ydata=ydata,
                maxfev=10000,
                p0=curr_params,
                bounds=(lb, ub),
                method='trf',
                full_output=True
            )
        except:
            print("mc =", _iter, ": failed")
            continue

        if ier in [1,2,3,4]:
            # a solution is found
            ypred= pap_simulation(xdata, *popt)
            curr_mse = np.sqrt(np.sum([(y-yhat)**2 for y,yhat in zip(ydata, ypred)])/len(xdata))
            mse.append(curr_mse)
            popt_list.append(popt)
            print("mc =", _iter, ": mse =", curr_mse)
        else:
            print("mc =", _iter, ": solution not found")
        _iter += 1

    if len(mse) == 0:
        raise RuntimeError("decomposition using one population model fails. check data.")
    else:
        mse_one = min(mse)
        popt_one = popt_list[mse.index(mse_one)]

    #----------------------------
    # run two population model
    #----------------------------
    print("\nrunning two population model...")
    use_one_population_model=False
    freq_s = 10**(uniform(np.log10(detection_limit), 0, size=mc_runs))
    deltag_s = lognormal(0, 1, size=mc_runs)
    deltag_r = lognormal(0, 1, size=mc_runs)
    hillc = uniform(1, 4, size=mc_runs)
    K50_s = exponential(10, size=mc_runs)
    K50_r = exponential(10, size=mc_runs)
    lb = [0.0,  0.0,  0.0,  1.0,  0.0,  0.0]
    ub = [1.0,  np.inf, np.inf, np.inf, np.inf, np.inf]
    mse = []
    popt_list = []
    _iter=1
    for p1, p2, p3, p4, p5, p6 in zip(freq_s, deltag_s, deltag_r, hillc, K50_s, K50_r):
        curr_params = [p1, p2, p3, p4, p5, p6]
        try:
            popt, pcov, infodict, mesg, ier = curve_fit(
                f=pap_simulation,
                xdata=xdata,
                ydata=ydata,
                maxfev=10000,
                p0=curr_params,
                bounds=(lb, ub),
                method='trf',
                full_output=True
            )
        except:
            print("mc =", _iter, ": failed")
            continue
        if ier in [1,2,3,4]:
            # a solution is found
            if popt[1]<popt[2]:
                tmp = popt[1]
                popt[1] = popt[2]
                popt[2] = tmp
            if popt[4]>popt[5]:
                tmp = popt[4]
                popt[4] = popt[5]
                popt[5] = tmp
            ypred= pap_simulation(xdata, *popt)
            curr_mse = np.sqrt(np.sum([(y-yhat)**2 for y,yhat in zip(ydata, ypred)])/len(xdata))
            mse.append(curr_mse)
            popt_list.append(popt)
            print("mc =", _iter, ": mse =", curr_mse)
        else:
            print("mc =", _iter, ": solution not found")
        _iter += 1

    if len(mse) == 0:
        raise RuntimeError("decomposition using two population model fails. check data.")
    else:
        mse_two = min(mse)
        popt_two = popt_list[mse.index(mse_two)]

    #-------------------------------------
    # compare one and two population model
    #-------------------------------------
    final_num_populations = 0
    if mse_one < mse_two + 0.001:
        # one population model
        results.append([strain, 'one', mse_one, 1.0, popt_one[0], np.NaN, popt_one[1], popt_one[2], np.NaN])
        use_one_population_model=True
        ypred = pap_simulation(xdata, *popt_one)
        final_num_populations = 1
    else:
        # two population model
        results.append([strain, 'two', mse_two] + list(popt_two))
        use_two_population_model=True
        ypred = pap_simulation(xdata, *popt_two)
        final_num_populations = 2
    assert final_num_populations in [1,2]
    print("\nmin mse_one =", mse_one, ", mser_two =", mse_two, ", %d population model is selected."%(final_num_populations))

    #-----
    # plot
    #-----
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6), sharex=True)

    # left panel:
    _ = ax[0].plot(xdata, ydata, 'o-', color='gray', markerfacecolor="none", markeredgecolor='gray', label='data')
    if final_num_populations == 1:
        _ = ax[0].plot(xdata, ypred, 'bx', label='fit', markersize=6)
    else:
        _ = ax[0].plot(xdata, ypred, 'bx', label='fit (freq R=%2.6f)'%(1-popt_two[0]))
    _ = ax[0].set_xlabel('Drug concentration')
    _ = ax[0].set_ylabel('log10(Survival)')
    _ = ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol=5)
    _ = ax[0].set_ylim([np.log10(detection_limit)-0.5,0.5])

    # right panel:
    if len(mic_to_plot) <=10:
        my_colors = sns.color_palette("tab10")
    elif len(mic_to_plot) <=20:
        my_colors = sns.color_palette("tab20")

    drug_concs_dense = np.linspace(0,np.max(drug_concs), 100)
    for k,mic in enumerate(mic_to_plot):
        if final_num_populations==1:
            maxg = popt_one[0]*(mic**popt_one[1]/(mic**popt_one[1] + popt_one[2]**popt_one[1]))
            if maxg <= maxg_ub:
                net_growth = maxg-popt_one[0]*(drug_concs_dense**popt_one[1]/(popt_one[2]**popt_one[1] + drug_concs_dense**popt_one[1]))
                _ = ax[1].plot(drug_concs_dense, net_growth, '-', color=my_colors[k], label="MIC=%2.2f"%mic)
        else:
            maxg = popt_two[1]*(mic**popt_two[3]/(mic**popt_two[3] + popt_two[4]**popt_two[3]))
            if maxg <= maxg_ub:
                net_growth_s = maxg-popt_two[1]*(drug_concs_dense**popt_two[3]/(popt_two[4]**popt_two[3]+drug_concs_dense**popt_two[3]))
                net_growth_r = maxg-popt_two[2]*(drug_concs_dense**popt_two[3]/(popt_two[5]**popt_two[3]+drug_concs_dense**popt_two[3]))
                if np.min(net_growth_r) > 0:
                    mic_r = np.inf # no MIC
                else:
                    mic_r = drug_concs_dense[np.argmin(np.abs(net_growth_r))],
                _ = ax[1].plot(drug_concs_dense, net_growth_s, '-', color=my_colors[k], label="MIC_s=%2.2f"%mic)
                _ = ax[1].plot(drug_concs_dense, net_growth_r, '--', color=my_colors[k], label="MIC_r=%2.2f"%mic_r)
    _ = ax[1].plot([0,np.max(drug_concs_dense)], [0,0], 'k--')
    _ = ax[1].set_ylabel('Net growth rate (1/hour)')
    _ = ax[1].set_xlabel('Drug concentration')
    _ = ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol=3)

    plt.tight_layout()
    plt.savefig("outputs/strain_%s.png" % (str(strain)), dpi=600)

df_results = pd.DataFrame(results, columns=['strain','pop_num','mse','freq_s','deltag_s','deltag_r','hillcoef','K50_s','K50_r'])
df_results.to_csv("decomposed_pap_parameters.csv")
