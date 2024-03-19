import matplotlib
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from nilearn import plotting

matplotlib.use('Agg')

sfc_nib = nib.load(snakemake.input.pscalar_sfc)
markers_nib = nib.load(snakemake.input.pscalar_markers)


sfc = sfc_nib.get_fdata().T 
markers = markers_nib.get_fdata().T

fig = plt.figure(figsize=(4,4))

print(sfc.shape)
print(sfc.max())
print(sfc.min())
print(np.isnan(sfc).sum())

print(markers.shape)
plotting.plot_markers(node_values=sfc,node_coords=markers,figure=fig)

fig.savefig(snakemake.output.png)

