wildcard_constraints:
    subject='[a-zA-Z0-9]+',
    app='snakebatch|freesurfer'


def get_zip_file(wildcards):
    #first check if this app is in the in_subj_zip list
    if wildcards.app in config['in_subj_zip']:
        subject = re.match('sub-([a-zA-Z0-9]+)',wildcards.file)[1] #parse subject id from the file
        return config['in_subj_zip'][wildcards.app].format(subject=subject)
    elif wildcards.app in config['in_agg_zip']:
        return config['in_agg_zip'][wildcards.app]
    else:
        print(f'ERROR: cannot find zip file for {wildcards.app} in config in_agg_zip or in_subj_zip')
        return None


rule get_from_zip:
    """ This is a generic rule to make any file within the {app} subfolder, 
        by unzipping it from a corresponding zip file"""
    input:
        zip=get_zip_file
    output:
        '{app}/{file}' # you could add temp() around this to extract on the fly and not store it
    shell:
        'unzip -d {wildcards.app} {input.zip} {wildcards.file}'

        
rule process_file:
    """ fake example of how you can make use of files that live inside the zipfile"""
    input: 
        pial='snakebatch/LOBE/derivatives/fmriprep/sub-{subject}/anat/sub-{subject}_acq-MPRvNavAvgEcho12_run-1_hemi-{hemi}_pial.surf.gii'
    output:
        postproc = bids(subject='{subject}',hemi='{hemi}',suffix='pial.surf.gii')
    shell: 
        'cp {input} {output}'
