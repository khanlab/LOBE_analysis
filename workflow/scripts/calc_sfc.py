import nibabel as nib
import numpy as np
from scipy.stats import spearmanr

struc_nib = nib.load(snakemake.input.pconn_struc)
func_nib = nib.load(snakemake.input.pconn_func)

#inputs are NxN
func_conn_matrix = func_nib.get_fdata()
struc_conn_matrix = struc_nib.get_fdata()

#create axes for output, we are using:
#  pscalar: ROW is scalars, COLUMN is parcels

# grab parcels axis from the input - can grab from axis=0 (row) or axis=1 (col)
parcel_axis = func_nib.header.get_axis(0)

# create new scalar axis, providing a name for each scalar
scalar_axis = nib.cifti2.ScalarAxis(name=['SFC']) 

N=func_conn_matrix.shape[0]

sfc = np.zeros((1,N))

#spearman rho over each row:
for i in range(N):
    func_row = func_conn_matrix[i,:]
    struc_row = struc_conn_matrix[i,:]

    #if all zeros in a row, then we set sfc to zero
    # since spearmanr will fail
    if func_row.sum() == 0 or struc_row.sum()==0:
        sfc[0,i] = 0
    else:
        sfc[0,i] = spearmanr(func_row,struc_row).statistic


#write out cifti
nib.Cifti2Image(sfc,header=(scalar_axis,parcel_axis)).to_filename(snakemake.output.pscalar_sfc)

