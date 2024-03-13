import pandas as pd
import numpy as np


#takes in confounds.tsv, adds spike regressors to, returns new confounds file
#read confounds table from fmriprep into a pandas dataframe
df = pd.read_table(snakemake.input.confounds_tsv)

fd_th = snakemake.params.spike_fd_th
dvars_th = snakemake.params.spike_dvars_th

outliers = (df.framewise_displacement > fd_th) | (df.std_dvars > dvars_th)

outliers = list(outliers[outliers].index)

n_volumes = len(df.index)

if outliers:
    spikes = np.zeros((n_volumes, len(outliers)))
    for i, outlier in enumerate(outliers):
        spikes[outlier, i] = 1.
        
    conf_spikes = pd.DataFrame(
        data=spikes, 
        columns=[f'motion_outlier_fd-{fd_th:01}_dvars-{dvars_th:01}_{i:02}' for i in range(len(outliers))]
        )

df_new = pd.concat((df,conf_spikes),axis=1)
n_spikes = len(outliers)


df_new.to_csv(snakemake.output.confounds_tsv,sep='\t')

