import pandas as pd
import nibabel as nib
from glob import glob
import numpy as np
import statsmodels.api as sp
#imported the full paths since it wasn't working otherwise - not sure why 
from statsmodels.stats.weightstats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path

df = pd.read_csv(snakemake.input.tsv,sep='\t')

group1_name=snakemake.wildcards.group1
group2_name=snakemake.wildcards.group2

group1 = df.query(snakemake.params.group1_query)
group2 = df.query(snakemake.params.group2_query)


group1['group']=group1_name
group2['group']=group2_name

combined_df = pd.concat([group1,group2])

Path(snakemake.output.plotsdir).mkdir(parents=True,exist_ok=True)

coltype = snakemake.params.coltype
print(df.columns)
columns = [col for col in df.columns if col.startswith(f'{coltype}_')]
print(columns)
for column_name in columns:
        
    x1 = group1[column_name]
    x2 = group2[column_name]
    
    tstats,pvals,dfs = ttest_ind(x1, x2, alternative='two-sided', usevar='pooled', weights=(None, None), value=0)
    plt.figure(figsize=(10,6))
    
    ax = sns.swarmplot(x='group',y=column_name,data=combined_df)
    
    pvals
    plt.title(f'{group1_name} vs {group2_name}, p-value={pvals}')

    out_png = f'{snakemake.output.plotsdir}/{column_name}.png'
    plt.savefig(out_png)
