def get_all_hippunfold_tsv(wildcards):
    tsvs={}
    for dataset in config['datasets']:
        tsvs[dataset]=expand(config['input_path']['hippunfold_vol_tsv'],
                        subject=subjects[dataset],
                        dataset=dataset)
    return tsvs


rule concat_hippunfold_vols:
    input:
        unpack(get_all_hippunfold_tsv)
    params:
        datasets=config['datasets']
    output:
        tsv = 'results/combined/hippunfold/group-all_desc-subfields_atlas-bigbrain_volumes.tsv'
    script:
        '../scripts/concat_tsv.py'
 
