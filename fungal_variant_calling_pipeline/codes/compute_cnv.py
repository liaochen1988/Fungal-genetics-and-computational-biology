import pandas as pd
import sys
import numpy as np
import os
from copy import deepcopy

df_cnv = None
current_folder = os.getcwd()
for _dir in os.listdir(current_folder+'/'+sys.argv[1]):
    # _dir is file name, not an absolute path
    if (len(sys.argv)==4 and not _dir.startswith(sys.argv[3])) or (not os.path.isdir(current_folder+'/'+sys.argv[1]+'/'+_dir)):
        continue

    # calculate reads coverage across the genome
    # genome_coverage.txt
    # Column0: chromosome (or entire genome)
    # Column1: depth of coverage from features in input file
    # Column2: number of bases on chromosome (or genome) with depth equal to column 1.
    # Column3: size of chromosome (or entire genome) in base pairs
    # Column4: fraction of bases on chromosome (or entire genome) with depth equal to column 1.
    mean_genome_coverage = None
    for _file in os.listdir(current_folder+'/'+sys.argv[2]+'/'+_dir):
        if _file.endswith(".genome_coverage.txt"):
            df_gc = pd.read_csv(current_folder+'/'+sys.argv[2]+'/'+_dir+'/'+_file, sep='\t', header=None, engine='python')
            df_gc = df_gc[df_gc[0]=='genome'] # use the statistics across the genome
            mean_genome_coverage = np.sum(df_gc[1]*df_gc[4])
            break
    if mean_genome_coverage is None:
        # could not find this file
        continue

    # calculate fature coverage on average
    # read_density.hist.txt
    # Column0: chromosome
    # Column1: start position
    # Column2: end position
    # Column3: gene name
    # Column4: depth of coverage
    # Column5: number of bases in the gene with depth equal to column 4.
    # Column6: length of gene
    # Column7: fraction of bases in the gene with depth equal to column 4.
    for _file in os.listdir(current_folder+'/'+sys.argv[2]+'/'+_dir):
        if _file.endswith(".read_density.hist.txt"):
            df_fc = pd.read_csv(current_folder+'/'+sys.argv[2]+'/'+_dir+'/'+_file, sep='\t', header=None, engine='python')
            df_fc = df_fc[df_fc[0] != 'all']
            df_fc[100] = df_fc[4]*df_fc[7]
            df_fc = df_fc[[0,3,100]]
            df_fc.columns = ['Chromosome','Feature','MeanCoverage']
            df_fc = df_fc.groupby(['Chromosome','Feature']).agg(np.sum).reset_index()
            df_fc['GenomeID'] = _file.rstrip('.read_density.hist.txt')
            df_fc['CopyNumber'] = df_fc['MeanCoverage']/mean_genome_coverage*2
            df_fc = df_fc[['GenomeID','Chromosome','Feature','CopyNumber']]
            break

    if df_cnv is None:
        df_cnv = deepcopy(df_fc)
    else:
        df_cnv = pd.concat([df_cnv, df_fc], ignore_index=True)

df_cnv.to_csv(current_folder+'/'+sys.argv[2]+'/ORF_v3_copy_number.csv', index=False)
