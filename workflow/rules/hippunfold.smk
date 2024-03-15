def get_all_hippunfold_tsv(wildcards):
    tsvs=[]
    for dataset in config['datasets']:
        tsvs.extend( 
            expand(config['input_path']['hippunfold_vol_tsv'],
                        subject=config['subjects'][dataset],
                        dataset=dataset))
    return tsvs


rule concat_hippunfold_vols:
    input:
        tsv_files = get_all_hippunfold_tsv
    output:
        concat_tsv = 'results/combined/hippunfold/group-all_desc-subfields_atlas-bigbrain_volumes.tsv'
    script:
        '../scripts/concat_tsv.py'
 
