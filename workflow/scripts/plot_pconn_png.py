import matplotlib
import matplotlib.pyplot as plt
import nibabel as nib
from nilearn import plotting

matplotlib.use('Agg')

cifti_nib = nib.load(snakemake.input.cifti_pconn)

conn_matrix = cifti_nib.get_fdata()

fig = plt.figure(figsize=(4,4))
plotting.plot_matrix(conn_matrix, figure=fig ,tri='lower', colorbar=False,
                         vmax=1, vmin=-1)

fig.savefig(snakemake.output.png)

