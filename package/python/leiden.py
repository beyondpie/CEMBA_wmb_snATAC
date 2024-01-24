# pyton version >= 3.10
import os
import sys
import math
import pickle
import random
from typing import Tuple, Union, Dict, List
import numpy as np
import itertools
from dataclasses import dataclass, field
import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
# https://matplotlib.org/stable/users/explain/interactive.html#controlling-interactive
# plt.ioff()
# colormaps
# https://matplotlib.org/stable/tutorials/colors/colormap-manipulation.html

from numba.core.errors import NumbaDeprecationWarning 
from numba.core.errors import NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

from snapatac2 import AnnData, AnnDataSet
from snapatac2._utils import is_anndata 
import snapatac2 as snap
from umap import UMAP



@dataclass
class LeidenSum:
    # clustering level
    cll:str
    # clustering id
    cid: str
    # barcodes
    barcodes: List[str]
    # resolutions
    rs: List[float]
    # silhouette scores
    sils: List[float]
    umap: np.ndarray
    name: str
    label: np.ndarray
    nsample: int = field(init = False)
    nc: int = field(init = False)
    def __post_init__(self):
        self.nsample = self.umap.shape[0]
        self.nc = len(self.rs)

# TODO: put it into one script for figure
@dataclass
class ScatterPlot:
    # enhance_title
    et: Dict[str, str] = field(default_factory = dict)
    ncol: int = 4
    width: int = 5
    nsample: int = 15000
    ncolors: int = 10
    palette: str = "viridis"
    # marker size
    ms: float = 2.0
    alpha: float = 1.0
    dpi: int = 200
    def __post_init__(self):
        self.et['fontweight'] = "bold"
        self.et["color"] = "red"

def cemba(
        adata: AnnData | np.ndarray,
        resolution: float,
        objective_function='modularity',
        n_repeat: int = 5,
        n_sample: int = 20000
) -> Tuple[float, float]:
    """CEMBA metrics for the selection of resolution of leiden algorithm.

    Detailes in "Cell clustering" of "Methods" in
    https://www.nature.com/articles/s41586-021-03604-1

    The metrics are all based on the so-called connectivity matrix M.
    - it's shape is (n_cell, n_cell)
    - M_ij \in (0, 1) is the fraction of cell i and cell j are in the same cluster
    during multiple round of repeats of leiden algorithms.

    In the original implementation, we use different random seeds for leiden
    algorithm as the multiple rounds of repeats. But we can also apply some random
    modifications of k-nearest neighbor (knn) graph, and then use leiden algorithm on the
    modified knn graph.

    Ideally, two cells are either in the same clusters or not. So most of M should be around
    zeros and ones. But if clustering is not that good, then M_{ij} may be around 0.5.

    Here we include three metrics:
    1. PAC, proportion of ambiguous clustering
    - it's \in (0,1)
    - the lower, the better
    - (\sum (M_{ij} < 0.95) - \sum (M_{ij} < 0.05) ) / n_cell ^2
    2. Dispersion coefficient (disp)
    - it's \in (0,1)
    - the higher, the better
    - \sum 4 * (M_{ij} - 0.5)^2 / n_cell ^2.

    3. Cumulative distribution function (cdf) curve 
    - draw the cdf based on the cdf of M.
    - since M_{ij} should be around 0 or 1 ideally. So the cdf curve would
        be have a relative flat region in the middle,
        and sharp increase around 0.0 and 1.0.

    During analysis, we use PAC and disp to choose the resolution, which
    will be the values of cemba_metrics.
    """
    membership = []
    for r in range(n_repeat):
        membership.append(snap.tl.leiden(
            adata,
            objective_function=objective_function, #pyright: ignore
            resolution=resolution,
            random_state=r, inplace=False
        ))
        partitions = np.array(membership).T
    return compute_metrics(
        partitions, #pyright: ignore
        nsample=n_sample, u1=0.05, u2=0.95)

def compute_metrics(
    partitions: np.ndarray,
    nsample: Union[int, None] = None,
    u1: float = 0.05,
    u2: float = 0.95
) -> Tuple[float, float]:
    """
    Parameters
    ----------
    partitions
        numpy.ndarray, dtype as np.unit, shape n_cell x n_times
    nsample
        int or None, used for downsampling cells, default is None.
    u1
        float, lower-bound for PAC, default is 0.05.
    u2
        float, upper-bound for PAC, default is 0.95.
    
    Returns
    -------
    Tuple[float, float]:
        (disp, PAC) in order.
    """
    
    # * check partitions
    ndim:int = partitions.ndim
    if ndim != 2:
        raise RuntimeError(
            f"partitions should have 2 instead of {ndim} dims.")
    # * to unsigned int64
    p = partitions.astype(np.int32)
    n_cell, n_times = p.shape
    # * downsample partitions if needed.
    if nsample and n_cell > nsample:
        print(f"Down sampling {n_cell} to {nsample}")
        random.seed(0)
        index = np.random.choice(n_cell, nsample, replace = False)
        p = p[index, :]
        n_cell, n_times = p.shape
    consensus = np.zeros((n_cell, n_cell), dtype = np.float16)
    for i in range(n_times):
        # print(f"{i+1} / {n_times} for consensus matrix")
        conn = cal_connectivity(p[:,i].tolist())
        consensus += conn
    consensus /= n_times
    disp: float = cal_dispersion(consensus)
    pac: float = cal_PAC(
        consensus, u1 = u1, u2 = u2)
    return (disp, pac)

def cal_connectivity(partition: list[int]) -> np.ndarray:
    """calculate connectivity matrix"""
    connectivity_mat = np.zeros((len(partition), len(partition)), dtype = bool)
    classN = max(partition)
    ## TODO: accelerate this
    for cls in range(int(classN + 1)):
        xidx = [i for i, x in enumerate(partition) if x == cls]
        iterables = [xidx, xidx]
        for t in itertools.product(*iterables):
            connectivity_mat[t[0], t[1]] = True
    """connectivity_mat = csr_matrix(connectivity_mat)"""
    return connectivity_mat

def cal_dispersion(consensus) -> float:
    """calculate dispersion coefficient

    Parameters
    ----------
    consensus
        Consensus matrix, shape (n_sample, n_sample). Each entry in the matrix is
        the fraction of times that two cells are clustered together.
    """
    n = consensus.shape[1]
    corr_disp = np.sum(
        4 * np.square(consensus - 0.5), dtype = np.float64) / (np.square(n))
    return corr_disp

def cal_PAC(consensus, u1, u2) -> float:
    """calculate PAC (proportion of ambiguous clustering)

    Parameters
    ----------
    consensus
        Consensus matrix, shape (n_sample, n_sample). Each entry in the matrix is
        the fraction of times that two cells are clustered together.
    """
    n = consensus.shape[0] ** 2
    PAC = ((consensus < u2).sum() - (consensus < u1).sum()) / n
    return PAC

def cal_silhouette(adata: AnnData | np.ndarray,
                   resolution: float,
                   membership: np.ndarray | None = None,
                   use_rep :str = "X_spectral",
                   rerun_leiden: bool = True,
                   leiden_key: str = "leiden",
                   idx: int = 0,
                   obj = 'modularity',
                   n_sample: int  = 20000,
                   random_state: int = 0,
                   logger = None) -> float:
    if membership is None:
        if not is_anndata(adata):
            raise RuntimeError("membership is None and adata is not AnnData.")
        if not rerun_leiden:
            membership = adata.obs.get(leiden_key, None) #pyright:ignore
            # we may have multiple clusters in leiden_key
            if membership and len(membership.shape) > 1:
                if membership.shape[1] == 1:
                    membership = membership[:, 0]
                else:
                    membership = membership[:, idx]
        if membership is None:
            membership = snap.tl.leiden(
                adata,
                objective_function= obj, #pyright:ignore
                resolution=resolution,
                random_state= random_state, inplace=False)
    if len(np.unique(membership)) < 2:
        if logger is not None:
            logger.warning(f"Under {resolution}, only one cluster found.")
        return -1.0
    from sklearn.metrics import silhouette_score
    data = adata.obsm[use_rep] if is_anndata(adata) else adata #pyright:ignore
    r = silhouette_score(
        X = data,
        labels = membership, #pyright: ignore
        sample_size = n_sample,
        random_state = random_state
    )
    return r #pyright: ignore

def init_LeidenSum_from_file(from_dir:str,
                             cll:str, cid: str,
                             prefix: str = "sa2_clustering") -> LeidenSum:
    nm = f"{prefix}_{cll}_{cid}"
    fnm = os.path.join(from_dir, f"{prefix}_{cll}_{cid}.pkl")
    print(f"Init LeidenSum from: {fnm}.")
    with open(fnm, 'rb') as f:
        id2info = pickle.load(f)
    t = LeidenSum(
        cll = cll,
        cid = cid,
        barcodes = id2info['barcode'],
        rs = id2info['leiden_r'],
        sils = id2info['sihouettes'],
        umap = id2info['umap'],
        # force to start from 1
        label = id2info['leiden'].astype(np.int64) + 1,
        name = nm
    )
    return t


def umap(
    adata: AnnData | AnnDataSet | np.ndarray,
    n_comps: int = 2,
    use_dims: int | list[int] | None = None,
    use_rep: str = "X_spectral",
    key_added: str = 'umap',
    random_state: int = 0,
    inplace: bool = True,
    **kwargs
) -> np.ndarray | None:
    """
    Parameters
    ----------
    adata
        The annotated data matrix.
    n_comps
        The number of dimensions of the embedding.
    use_dims
        Use these dimensions in `use_rep`.
    use_rep
        Use the indicated representation in `.obsm`.
    key_added
        `adata.obs` key under which to add the cluster labels.
    random_state
        Random seed.
    inplace
        Whether to store the result in the anndata object.

    Returns
    -------
    np.ndarray | None
        if `inplace=True` it stores UMAP embedding in
        `adata.obsm["X_`key_added`"]`.
        Otherwise, it returns the result as a numpy array.
    """
    data = adata.obsm[use_rep] if is_anndata(adata) else adata

    if use_dims is not None:
        data = data[:, :use_dims] if isinstance(use_dims, int) else data[:, use_dims]

    umap = UMAP(
        random_state=random_state,
        n_components=n_comps,
        **kwargs
    ).fit_transform(data)
    if inplace:
        adata.obsm["X_" + key_added] = umap
    else:
        return umap


def draw_umap(t: LeidenSum,
              # scatter formatter
              scf: ScatterPlot,
              colors: List[str] = ["black"],
              outdir: str|None = None) -> Figure:
    ncol:int = scf.ncol
    nrow:int = math.ceil(t.nc / ncol)
    fig, axs = plt.subplots(
        nrow, ncol,
        figsize = (scf.width*nrow, scf.width*ncol),
        constrained_layout = True,
        squeeze = False)
    which_max = np.argmax(t.sils)
    # prepare data
    index = range(t.nsample)
    if t.nsample > scf.nsample:
        print(f"Down sample data to size of {scf.nsample} from {t.nsample}.")
        index = np.random.choice(range(t.nsample),
                                 size = scf.nsample, replace = False)
    x = t.umap[index, 0]
    y = t.umap[index, 1]
    for i in range(nrow):
        for j in range(ncol):
            ax = axs[i, j]
            ith = i*ncol + j
            # get raw nclusters
            nc = len(np.unique(t.label[:,ith]))
            label = t.label[index,ith].astype(np.int64)
            ulabel_ = np.unique(label)

            # set colors
            if len(colors) >= len(ulabel_):
                colormap = ListedColormap(colors)
            else:
                colormap = matplotlib.colormaps[
                    scf.palette].resampled(len(ulabel_))
            sc = ax.scatter(x = x, y = y,
                            c = label,
                            s = scf.ms,
                            alpha = scf.alpha,
                            edgecolors = 'none',
                            cmap = colormap)
            title = f"r_{t.rs[ith]}, sil_{round(t.sils[ith],2)}, nc_{nc}"
            ax.set_title(title)
            if ith == which_max:
                ax.set_title(title, fontdict = scf.et)
            ax.legend(handles = sc.legend_elements()[0],
                      labels = ulabel_.tolist())
    print("Generated all the umaps.")
    if outdir is not None:
        print(f"Output to: {outdir}.")
        os.makedirs(outdir, exist_ok = True)
        outf = os.path.join(outdir, f"sa2_sum_{t.cll}_{t.cid}_umap.png")
        fig.savefig(fname = outf, dpi = scf.dpi, bbox_inches = 'tight')
    return fig
