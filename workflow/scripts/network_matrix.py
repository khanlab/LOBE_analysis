#get atlas
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from nilearn import plotting


import matplotlib
matplotlib.use('Agg')


#load schaefer atlas networks

#'resources/schaefer_2018/Schaefer2018_300Parcels_7Networks_order.txt'

df_networks = pd.read_table(snakemake.input.network_txt,
                                sep='\t',header=None,
                                names=['lbl','network','x','y','z','_'])


networks7 = list()
for network in df_networks.network:
    hemi = network.split('_')[1]
    net = network.split('_')[2]
    
    networks7.append(f'{hemi}_{net}')
    
df_networks['net7']=networks7

net7_names = df_networks['net7'].unique()


# now get network connectivity matrix by averaging over network rois

# load the Schaefer connectivity first
conn = np.loadtxt(snakemake.input.conn_txt)
n_rois = conn.shape[0]
#print(conn)

net_conn = np.zeros((len(net7_names),len(net7_names)))

for i,net_i in enumerate(net7_names):
    for j,net_j in enumerate(net7_names):

        i_mask = df_networks['net7'] == net_i
        j_mask = df_networks['net7'] == net_j

        mask_rows = np.zeros(conn.shape, dtype='bool')
        mask_cols = np.zeros(conn.shape,dtype='bool')

        mask_rows[i_mask,:] = True
        mask_cols[:,j_mask] = True

        mask = mask_rows & mask_cols

        net_conn[i,j] = conn[mask].mean()
        
         

from nilearn import plotting

fig = plt.figure(figsize=(4,4))
plotting.plot_matrix(net_conn, figure=fig ,tri='diag', colorbar=False,
                         vmax=1, vmin=-1,labels=net7_names)


fig.savefig(snakemake.output.png)
np.savetxt(snakemake.output.txt,net_conn)

