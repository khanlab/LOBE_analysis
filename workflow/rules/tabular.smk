#rule import_hcp_tabular:
#    input:
#        ndar_txt=lambda wildcards: Path(in_root) / config["in_tabular"][wildcards.dataset],
#    params:
#        subjects=lambda wildcards: subjects[wildcards.dataset]
#    output:
#        tsv="resources/dataset-{dataset,HCP}_tabular.tsv",
#    script:
#        "../scripts/import_hcp_tabular.py"

rule import_lobe_tabular:
    input:
        tsv=lambda wildcards: Path(in_root) / config["in_tabular"][wildcards.dataset],
    params:
        subjects=lambda wildcards: subjects[wildcards.dataset]
    output:
        tsv="resources/dataset-{dataset,LOBE}_tabular.tsv",
    script:
        "../scripts/import_lobe_tabular.py"

rule merge_tabular:
    input:
        subjects_tsvs=expand("resources/dataset-{dataset}_subjects.tsv",
            dataset=config['datasets']),
        tabular_tsvs=expand("resources/dataset-{dataset}_tabular.tsv",
            dataset=config['datasets'])
    params:
        datasets = config['datasets'],
        subjects=subjects
    output:
        tsv="resources/merged_tabular.tsv"
    script: "../scripts/merge_tabular.py"

rule add_conn_to_tabular:
    input:	
        get_dwi_targets(),
        get_func_targets(),
        get_sfc_targets(),
        tsv="resources/merged_tabular.tsv",
    params:
        subjects=subjects
    output:
        tsv="resources/merged_tabular_withconn.tsv"
    notebook:
        '../notebooks/add_conn_to_tabular.py.ipynb'	
