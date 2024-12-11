import numpy as np
import scanpy as sc


# Preprocessing
# filter out genes only expressed in a few spots
# filter out spots where only a few genes are detected
def spark_filter(adata, min_spot_percentage=0.1, min_gene_counts=10, inplace=True):

    if min_spot_percentage > 0:
        sc.pp.filter_genes(adata, min_cells=np.floor(min_spot_percentage * adata.shape[0]), inplace=inplace)
    if min_gene_counts > 0:
        sc.pp.filter_cells(adata, min_counts=min_gene_counts, inplace=inplace)
    if not inplace:
        return adata


# Calculate the parameters for gaussian kernel and cosine kernel
def cal_kernel_params(mtx=None, dist_mtx=None, min_idx=3, counts=5):

    if dist_mtx is not None:
        dmin = np.min(dist_mtx[dist_mtx >= 1e-8]) / 2
        dmax = np.max(dist_mtx[dist_mtx >= 1e-8]) * 2
    else:
        mtx = np.array(mtx)
        diff = mtx[:, np.newaxis, :] - mtx[np.newaxis, :, :]
        dist_mtx = np.sqrt(np.sum(diff ** 2, axis=2))

        dmin = np.min(dist_mtx[dist_mtx >= 1e-8]) / 2
        dmax = np.max(dist_mtx[dist_mtx >= 1e-8]) * 2
    params = 10 ** np.linspace(np.log10(dmin), np.log10(dmax), 11)[min_idx: (min_idx + counts)]

    return list(params)


# Gaussian kernel
def gaussian_kernel(mtx=None, dist_mtx=None, sigma=0.1):

    if dist_mtx is not None:
        dist_square = dist_mtx ** 2
    else:
        mtx = np.array(mtx)
        diff = mtx[:, np.newaxis, :] - mtx[np.newaxis, :, :]
        dist_square = np.sum(diff ** 2, axis=2)
    kernel_mtx = np.exp(-dist_square / (2 * sigma ** 2))

    return kernel_mtx


# Cosine kernel
def cosine_kernel(mtx=None, dist_mtx=None, phi=0.1):

    if dist_mtx is None:
        mtx = np.array(mtx)
        diff = mtx[:, np.newaxis, :] - mtx[np.newaxis, :, :]
        dist_mtx = np.sqrt(np.sum(diff ** 2, axis=2))
    kernel_mtx = np.cos(2 * np.pi * dist_mtx / phi)

    return kernel_mtx





