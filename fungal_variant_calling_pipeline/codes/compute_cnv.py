import pandas as pd
import numpy as np
import os
from copy import deepcopy

df_cnv = None
for _dir in os.listdir(sys.argv[1]):

    # calculate genome coverage on average
    mean_genome_coverage = None
    for file in os.listdir(_dir):
        if file.endswith(".genome_coverage.txt"):
            df_gc = pd.read_csv('%s/%s'%(_dir, file), sep='\t', header=None, engine='python')
            df_gc = df_gc[df_gc[0]=='genome'] # remove mitochondrial
            # columns: chromosome, coverage, counts of bases that have coverage in the second column, genome size, percentage of genome that has the specific coverage depth
            mean_genome_coverage = np.sum(df_gc[1]*df_gc[4])
            break
    if mean_genome_coverage is None:
        # could not find this file
        continue
        
    # calculate fature coverage on average
    for file in os.listdir(_dir):
        if file.endswith(".read_density.hist.txt"):
            df_fc = pd.read_csv('%s/%s'%(_dir, file), sep='\t', header=None, engine='python')
            df_fc = df_fc[df_fc[0] != 'all']
            # columns: chromosome, start pos, end pos, feature name, coverage, counts of bases in this feature that have the coverage in the previous column, length of the feature, percentage of the feature that has the specific coverage depth
            df_fc[100] = df_fc[4]*df_fc[7]
            df_fc = df_fc[[0,3,100]]
            df_fc.columns = ['Chromosome','Feature','MeanCoverage']
            df_fc = df_fc.groupby(['Chromosome','Feature']).agg(np.sum).reset_index()
            df_fc['GenomeID'] = file.split('_')[0]
            df_fc['CopyNumber'] = df_fc['MeanCoverage']/mean_genome_coverage*2
            break

    if df_cnv is None:
        df_cnv = deepcopy(df_fc)
    else:
        df_cnv = pd.concat([df_cnv, df_fc], ignore_index=True)
        
df_cnv.to_csv('%s/orf_copy_number.csv' % (sys.argv[2]), index=False)