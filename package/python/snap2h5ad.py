"""
Save mat from snap file to h5ad file.
"""
import anndata as ad
import h5py
import numpy as np
from scipy import sparse

def snapbmat2h5ad(snap_file, bin_size:int = 5000,
                  prefix:str = '',
                  outfnm:str = '') -> ad.AnnData:
    """Save bmat from snap file to h5ad.
    """
    with h5py.File(snap_file, mode = 'r') as f:
        barcode = [b.decode('utf-8') for b in f['/BD/name']]
        chrom = [i.decode('utf-8')
                     for i in f[f'/AM/{bin_size}/binChrom']]
        start = f[f'/AM/{bin_size}/binStart'][:]
        end = start + bin_size - 1
        name = [f"{i}:{s}-{t}" for i, s, t in
                    zip(chrom, start, end)]
        idx = np.array(f[f'/AM/{bin_size}/idx'],
                            dtype = np.uintc)
        idy = np.array(f[f'/AM/{bin_size}/idy'],
                            dtype = np.uintc)
        count = np.array(f[f'/AM/{bin_size}/count'],
                              dtype = np.float32)
    cscmat = sparse.csc_matrix(
        (count, (idx - 1, idy - 1)),
        shape = (len(barcode), len(name)), dtype = np.float32)
    adata = ad.AnnData(X = cscmat)
    if len(prefix) > 0:
        fullname = [f"{prefix}.{b}" for b in barcode]
    else:
        fullname = barcode
    adata.obs_names = fullname
    adata.var_names = name
    if len(outfnm) > 0 :
        adata.write(outfnm, compression = "gzip")
    return adata


def snapgmat2h5ad(snap_file, prefix:str = '',
                  outfnm:str = '') -> ad.AnnData:
    """Save gmat from snap file to h5ad.
    """
    with h5py.File(snap_file, mode = 'r') as f:
        barcode = [b.decode('utf-8') for b in f['/BD/name']]
        genm = [g.decode('utf-8') for g in f['/GM/name']]
        idx = np.array(f["/GM/idx"], dtype = np.uintc)
        idy = np.array(f["/GM/idy"], dtype = np.uintc)
        count = np.array(f["/GM/count"], dtype = np.float32)
    cscmat = sparse.csc_matrix(
        (count, (idx -1, idy-1)),
        shape = (len(barcode), len(genm)),
        dtype = np.float32
    )
    adata = ad.AnnData(X = cscmat)
    if len(prefix) > 0:
        fullname = [f"{prefix}.{b}" for b in barcode]
    else:
        fullname = barcode
    adata.obs_names = fullname
    adata.var_names = genm
    if len(outfnm) > 0 :
        adata.write(outfnm, compression = "gzip")
    return adata
