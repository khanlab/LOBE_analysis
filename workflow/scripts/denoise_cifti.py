import json
import nibabel as nib
import pandas as pd
import numpy as np
from snakebids import write_derivative_json
from nilearn.image import clean_img 

#read confounds table from fmriprep into a pandas dataframe
df = pd.read_table(snakemake.input.confounds_tsv)

#get TR (needed for filtering):     
with open(snakemake.input.json) as f:
    t_r = json.load(f)['RepetitionTime']
 
#load 4d bold
cifti_nib = nib.load(snakemake.input.cifti)

bold_data = cifti_nib.get_fdata() # time x voxels
(n_timesteps,n_voxels) = bold_data.shape

#transpose and reshape to be like 4d image
bold_nib = nib.Nifti1Image(bold_data.T.reshape((n_voxels,1,1,n_timesteps)),affine=np.eye(4))

#get denoise parameters:
denoise_dict = snakemake.params.denoise_params

#get specified confounds from the dataframe
df_list = list()
df_list.append( df.filter(items = denoise_dict['confounds_name']))

for term in denoise_dict['confounds_like']:
    df_list.append( df.filter(like = term))

df_filtered = pd.concat(df_list,axis=1)

print(df_filtered.columns)

confounds = df_filtered.to_numpy()

#print('before removing non-finite')
#print(confounds)
#replace non-finite values with zeros in the confounds
confounds[np.logical_not(np.isfinite(confounds))] = 0


#use nilearn to clean
cleaned = clean_img(imgs=bold_nib,
                    confounds=confounds,
                    t_r=t_r,
                    **denoise_dict['clean_img_opts'])

#save to nifti
cleaned_data = cleaned.get_fdata().reshape((n_voxels,n_timesteps)).T
cifti_nib_out = nib.Cifti2Image(cleaned_data,header=cifti_nib.header)

cifti_nib_out.to_filename(snakemake.output.cifti)


