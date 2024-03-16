rule reslice_mask:
    input: 
        boldmask = lambda wildcards: config['input_path']['bold_mask'][wildcards.dataset],
        t1mask = lambda wildcards: config['input_path']['t1_mask'][wildcards.dataset],
    output:
        mask = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',suffix='brain_mask.nii.gz')
    shell:
        'c3d {input.boldmask} {input.t1mask} -interpolation NearestNeighbor  -reslice-identity -o {output}'

rule denoise:
    input: 
        nii = lambda wildcards: config['input_path']['bold_nii'][wildcards.dataset],
        json = lambda wildcards: config['input_path']['bold_json'][wildcards.dataset],
        confounds_tsv = lambda wildcards: config['input_path']['bold_confounds'][wildcards.dataset],
	    mask_nii = rules.reslice_mask.output.mask
    params:
        denoise_params = lambda wildcards: config['fmri']['denoise'][wildcards.denoise],
    output: 
        nii = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',suffix='bold.nii.gz'),
        json = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',suffix='bold.json')
    group: 'subj'
    threads: 8
    resources:
        mem_mb='32000'
    script: '../scripts/denoise.py'

  
rule smooth:
    input:
        nii = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',suffix='bold.nii.gz'),
        json = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',suffix='bold.json')
    params:
        fwhm = lambda wildcards: float(wildcards.fwhm)
    output:
        nii = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',suffix='bold.nii.gz'),
        json = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',suffix='bold.json')
    group: 'subj'
    script: '../scripts/smooth.py'

 
rule schaefer_connectivity:
    input:
        nii = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',suffix='bold.nii.gz'),
    params:
        n_rois = 300,
        yeo_networks = 7,
        data_dir = 'resources',
        conn_measure = 'correlation'
    group: 'subj'
    output:
        txt = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',atlas='{atlas,schaefer}',suffix='conn.txt'),
        png = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',atlas='{atlas,schaefer}',suffix='conn.png'),
    script: '../scripts/connectivity_matrix.py'


