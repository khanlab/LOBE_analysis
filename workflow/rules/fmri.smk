
           
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

 
rule schaefer_network:
    """ Uses the Schaefer connectivity and 7-network labels to construct 14x14 network connectivity matrices, by averaging values across the network parcels """
    input:
        conn_txt = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',atlas='{atlas}',suffix='conn.txt'),
        network_txt = 'resources/schaefer_2018/Schaefer2018_300Parcels_7Networks_order.txt' #this should be downloaded when the prev rule is first run
    group: 'subj'
    output:
        txt = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',atlas='{atlas,schaefer}',network='Yeo7',suffix='conn.txt'),
        png = bids(root=root,subject='{subject}',task='{task}',denoise='{denoise}',space='{space}',fwhm='{fwhm}',atlas='{atlas,schaefer}',network='Yeo7',suffix='conn.png'),
    script: '../scripts/network_matrix.py'
