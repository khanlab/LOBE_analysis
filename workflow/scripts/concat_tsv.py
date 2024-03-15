import pandas as pd

df_list=[]
for dataset in snakemake.params.datasets:

    df = pd.concat([pd.read_table(in_tsv) for in_tsv in snakemake.input[dataset]])
    df['dataset'] = dataset
    df_list.append(df)


pd.concat(df_list).to_csv(snakemake.output.tsv, sep="\t", index=False)
