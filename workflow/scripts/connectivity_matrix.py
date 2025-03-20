#get atlas

import matplotlib.pyplot as plt
import numpy as np
from nilearn import plotting
from nilearn import datasets
from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure


import matplotlib
matplotlib.use('Agg')


correlation_measure = ConnectivityMeasure(kind=snakemake.params.conn_measure)


#schaefer
schaefer = datasets.fetch_atlas_schaefer_2018(n_rois=snakemake.params.n_rois,
                                                yeo_networks=snakemake.params.yeo_networks,
                                                data_dir=snakemake.params.data_dir)
atlas_maps = schaefer.maps
atlas_labels = schaefer.labels

masker = NiftiLabelsMasker(labels_img=atlas_maps,verbose=5)

timeseries = masker.fit_transform(snakemake.input.nii)
print(timeseries.shape)
correlation_matrix = correlation_measure.fit_transform([timeseries])[0]

fig = plt.figure(figsize=(4,4))
plotting.plot_matrix(correlation_matrix, figure=fig ,tri='lower', colorbar=False,
                         vmax=1, vmin=-1)

fig.savefig(snakemake.output.png)
np.savetxt(snakemake.output.txt,correlation_matrix)


