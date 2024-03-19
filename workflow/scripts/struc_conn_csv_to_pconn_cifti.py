import numpy as np
import nibabel as nib


ref_cifti_nib = nib.load(snakemake.input.ref_cifti_pconn)

conn_matrix = np.loadtxt(snakemake.input.conn_csv,delimiter=',')

nib.Cifti2Image(conn_matrix,header=ref_cifti_nib.header).to_filename(snakemake.output.cifti_pconn)

