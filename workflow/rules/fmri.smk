
           
wildcard_constraints: 
    subject = '[a-zA-Z0-9]+',
    task = '[a-zA-Z0-9]+',
    denoise = '[a-zA-Z0-9]+',
    space = '[a-zA-Z0-9]+',
    atlas = '[a-zA-Z0-9]+',
    fwhm = '[a-zA-Z0-9]+'

         
root = os.path.join(config["root"], "{dataset}","fmri")

def get_fmri_targets():
    targets=[]
    for dataset in config['datasets']:
        targets.extend( 
            expand(bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',
                space='{space}',fwhm='{fwhm}',atlas='{atlas}',network='Yeo7',suffix='conn.png'),
        subject=config['subjects'][dataset],
        dataset=dataset,
        task=config['fmri']['task'],
        denoise=config['fmri']['denoise'].keys(),
        space=config['fmri']['space'],
        fwhm=config['fmri']['fwhm'],
        atlas=config['fmri']['atlas']
        ))
    return targets

      


rule all_fmri:
    input:
        get_fmri_targets()


rule map_bold_to_surface_fsLR:
    input:
        bold_preproc = lambda wildcards: config['input_path']['bold_nii'][wildcards.dataset],
        mid_surf=lambda wildcards: config['input_path']['surf_gii_mni'].format(surf='midthickness',**wildcards),
        white_surf=lambda wildcards: config['input_path']['surf_gii_mni'].format(surf='white',**wildcards),
        pial_surf=lambda wildcards: config['input_path']['surf_gii_mni'].format(surf='pial',**wildcards),
    output:
        metric=bids(root=root,datatype='func',hemi='{hemi}',desc='preproc',space='{space}',den='32k',task='{task}',suffix='bold.dtseries.func.gii',
                **config['subj_wildcards'])
    shell:
        'wb_command -volume-to-surface-mapping {input.bold_preproc} {input.mid_surf}  {output.metric} -ribbon-constrained {input.white_surf} {input.pial_surf}'


rule create_bold_cifti:
    input:
        left_metric=bids(root=root,datatype='func',hemi='L',desc='preproc',space='{space}',den='32k',task='{task}',suffix='bold.dtseries.func.gii',
                **config['subj_wildcards']),
        right_metric=bids(root=root,datatype='func',hemi='R',desc='preproc',space='{space}',den='32k',task='{task}',suffix='bold.dtseries.func.gii',
                **config['subj_wildcards']),
    output:
        cifti=bids(root=root,datatype='func',desc='preproc',space='{space}',den='32k',task='{task}',suffix='bold.dtseries.nii',
                **config['subj_wildcards']),
    shell:
        'wb_command -cifti-create-dense-timeseries {output.cifti} -left-metric {input.left_metric} -right-metric {input.right_metric} '

rule denoise_cifti:
    input: 
        cifti=rules.create_bold_cifti.output.cifti,
        json = lambda wildcards: config['input_path']['bold_json'][wildcards.dataset],
        confounds_tsv = lambda wildcards: config['input_path']['bold_confounds'][wildcards.dataset],
    params:
        denoise_params = lambda wildcards: config['fmri']['denoise'][wildcards.denoise],
    output: 
        cifti=bids(root=root,datatype='func',desc='preproc',space='{space}',den='32k',task='{task}',denoise='{denoise}',suffix='bold.dtseries.nii',
                **config['subj_wildcards']),
    group: 'subj'
    threads: 8
    resources:
        mem_mb='32000'
    script: '../scripts/denoise_cifti.py'


rule smooth_cifti:
    input:
        cifti=rules.denoise_cifti.output.cifti,
        left_surf=lambda wildcards: config['input_path']['surf_gii_t1'].format(surf='midthickness',hemi='L',**wildcards),
        right_surf=lambda wildcards: config['input_path']['surf_gii_t1'].format(surf='midthickness',hemi='R',**wildcards),
    params:
        fwhm = lambda wildcards: float(wildcards.fwhm)
    output:
        cifti=bids(root=root,datatype='func',desc='preproc',space='{space}',den='32k',task='{task}',denoise='{denoise}',fwhm='{fwhm}',suffix='bold.dtseries.nii',
                **config['subj_wildcards']),
    shell:
        'wb_command -cifti-smoothing {input.cifti} {params.fwhm} {params.fwhm} COLUMN '
        ' {output.cifti} -fwhm -left-surface {input.left_surf} -right-surface {input.right_surf}' 


rule parcellate_bold:
    input: 
        cifti_dtseries=rules.smooth_cifti.output.cifti,
        cifti_dlabel=lambda wildcards: config['atlas'][wildcards.atlas]
    params:
        exclude_opt='-exclude-outliers {nstdev} {nstdev}'.format(nstdev=config['fmri']['parcellation']['n_stdevs_exclude']) if config['fmri']['parcellation']['do_exclude_outliers'] else ''
    output:
        cifti_ptseries=bids(root=root,datatype='func',desc='preproc',space='{space}',den='32k',task='{task}',denoise='{denoise}',fwhm='{fwhm}',atlas='{atlas}',suffix='bold.ptseries.nii',
                **config['subj_wildcards']),
    shell:
        'wb_command -cifti-parcellate {input.cifti_dtseries} {input.cifti_dlabel} '
        ' COLUMN {output.cifti_ptseries} {params.exclude_opt}'
       
rule correlate_parcels:
    input: 
        cifti=rules.parcellate_bold.output.cifti_ptseries
    params:
        fisher_z='-fisher-z' if config['fmri']['correlation']['do_fisher_z'] else ''
    output:
        cifti=bids(root=root,datatype='func',desc='preproc',space='{space}',den='32k',task='{task}',denoise='{denoise}',fwhm='{fwhm}',atlas='{atlas}',suffix='bold.pconn.nii',
                **config['subj_wildcards']),
    shell:
        'wb_command -cifti-correlation {input.cifti} {output.cifti} {params.fisher_z} ' 


rule plot_pconn_png:
    """generic rule for plotting pconn cifti files"""
    input:
        cifti_pconn='{prefix}.pconn.nii'
    output:
        png='{prefix}.pconn.png'
    script:
        '../scripts/plot_pconn_png.py'

 
rule struc_conn_csv_to_pconn_cifti:
    input:
        #for reference pconn, pick the first task, denoise, fwhm values from config (any would do)
        ref_cifti_pconn=bids(root=root,datatype='func',desc='preproc',space='{space}',den='32k',task=config['fmri']['task'][0],denoise=next(iter(config['fmri']['denoise'])),fwhm=config['fmri']['fwhm'][0],atlas='{atlas}',suffix='bold.pconn.nii',
                **config['subj_wildcards']),
        conn_csv=bids(
            root=mrt_root,
            datatype='tractography',
            atlas='{atlas}',
            suffix='conn.csv',
            **config['subj_wildcards'],
        )
    output:
        cifti_pconn=bids(
            root=mrt_root,
            datatype='tractography',
            space='{space}',
            den='32k',
            atlas='{atlas}',
            suffix='struc.pconn.nii',
            **config['subj_wildcards'],
        )
    script:
        '../scripts/struc_conn_csv_to_pconn_cifti.py'



    

