"""
After running nmf independently, we read the results and summize them
based on Yang's nmf script.
"""

import argparse

parser = argparse.ArgumentParser(description="Run NMF using sklearn.")
parser.add_argument("-i", "--inputPrefix", type=str, dest="input_prefix")
parser.add_argument("-r", "--rank", type=int, dest="rank")
parser.add_argument("-d", "--defaultRandomState", type = int,
                    dest = "default_random_state")
parser.add_argument("-s", "--startRandomState", type = int,
                    dest = "start_random_state")
parser.add_argument("-e", "--endRandomState", type = int,
                    dest = "end_random_state")
parser.add_argument("-o", "--outPrefix", type = str,
                    dest="out_prefix")

args = parser.parse_args()

from warnings import warn
from os.path import exists, basename
import numpy as np
from scipy import io
from scipy.cluster.hierarchy import linkage, leaves_list, cophenet
from scipy.spatial.distance import squareform
from scipy.sparse import save_npz, load_npz
from sklearn.decomposition import NMF
import fastcluster as fc
from time import perf_counter as pc
import matplotlib.pyplot as plt

plt.switch_backend("agg")
try:
    from matplotlib.pyplot import imshow, set_cmap
except ImportError as exc:
    warn("Matplotlib must be installed.")

import itertools
import math

def run():
    """Run post analysis to summaize nmf."""
    input_prefix = args.input_prefix
    out_prefix = args.out_prefix
    rank = args.rank
    default_random_state = args.default_random_state
    start_random_state = args.start_random_state
    end_random_state = args.end_random_state
    start_time = pc()
    post_nmf_2(input_prefix = input_prefix,
               n_mods = rank,
               out_prefix = out_prefix,
               start_random_state = start_random_state,
               end_random_state = end_random_state,
               default_random_state = default_random_state)
    end_time = pc()
    print("Used (secs): ", end_time - start_time)


def read_npz(inputf):
    """
    Read snATAC data in npz format. The matrix's shape is ### (bins) x ### (samples).
    It contains only positive data (boolean).

    Return the dense data matrix.
    """
    V = load_npz(inputf)
    V = V.tocsr()
    return V


def read_mtx(inputf):
    print("=== read in MM format files %s ===" % inputf)
    data = io.mmread(inputf)
    data = data.tocsr()
    return data


def save_npz_gi(V, xgi, ygi, prefix):
    print("=== save npz format with xgi and ygi===")
    npz_file = ".".join([prefix, "npz"])
    save_npz(npz_file, V)
    xgi_file = ".".join([prefix, "xgi"])
    np.savetxt(xgi_file, xgi, fmt="%s", delimiter="\t")
    ygi_file = ".".join([prefix, "ygi"])
    np.savetxt(ygi_file, ygi, fmt="%s", delimiter="\t")


def read_xgi(xgi):
    fname = xgi
    with open(fname) as f:
        xgi_content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
    xgi_content = [x.strip() for x in xgi_content]
    return xgi_content


def read_ygi(ygi):
    fname = ygi
    with open(fname) as f:
        ygi_content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
    ygi_content = [x.strip() for x in ygi_content]
    return ygi_content


def write_coo(coo, prefix):
    coo = coo.tocoo()
    coo_file = ".".join([prefix, "mtx"])
    f = open(coo_file, "w")
    for row, col, value in zip(coo.row, coo.col, coo.data):
        out = "{0}\t{1}\t{2}\n".format(row, col, value)
        f.write(out)


def run_one(V, rank, outPrefix):
    print("=== perform NMF on rank %d ===" % rank)
    model = NMF(n_components=rank, init="random", random_state=0, verbose=True)
    W = model.fit_transform(V)
    H = model.components_
    saveW(outPrefix, W)
    saveH(outPrefix, H)
    print(rank)
    print(model.reconstruction_err_)
    print(model.n_iter_)


def saveC(prefix, X):
    print("=== write matrix C ===")
    fileN = [prefix, "C", "mx"]
    fileN = ".".join(fileN)
    np.savetxt(fileN, X, fmt="%g", delimiter="\t")


def saveH(prefix, X):
    print("=== write matrix H ===")
    fileN = [prefix, "H", "mx"]
    fileN = ".".join(fileN)
    np.savetxt(fileN, X, fmt="%g", delimiter="\t")


def saveW(prefix, X):
    print("=== write matrix W ===")
    fileN = [prefix, "W", "mx"]
    fileN = ".".join(fileN)
    np.savetxt(fileN, X, fmt="%g", delimiter="\t")


def plotH(prefix, X):
    """
    Plot reordered consensus matrix.

    :param C: Reordered consensus matrix.
    :type C: numpy.ndarray`
    :param rank: Factorization rank.
    :type rank: `int`
    """
    fig = plt.figure(figsize=(5, 5), dpi=100)
    imshow(X)
    set_cmap("RdBu_r")
    fileN = [prefix, "H", "png"]
    fileN = ".".join(fileN)
    fig.savefig(fileN)


def cal_sparseness(X):
    print("=== calculate sparseness ===")
    vec = list(np.concatenate(X))
    absSum = np.sum(np.abs(vec))
    n = np.prod(X.shape)
    squareSum = np.sum(np.square(vec))
    numerator = np.sqrt(n) - (absSum / np.sqrt(squareSum))
    denominator = np.sqrt(n) - 1
    sparseness = numerator / denominator
    return sparseness


def cal_rss_mse(W, H, V):
    """Residual Sum of Squares (RSS) & Mean Square Error (MSE)"""
    print("=== calculate Residual Sum of Squares (RSS) & Mean Square Error (MSE) ===")
    residualSquare = np.square(W.dot(H) - V)
    rss = np.sum(residualSquare)
    mse = np.mean(residualSquare)
    out = [rss, mse]
    return out


def cal_mse(W, H, V):
    """# The mean square error"""
    print("=== calculate Mean Square Error (MSE) ===")
    mse = np.mean(np.square(W.dot(H) - V))
    return mse


def cal_evar(rss, V):
    print("=== calculate evar ===")
    evar = 1 - (rss / np.sum(V.data ** 2))
    return evar


def cal_featureScore_kim(W):
    """extract feature from W"""
    print("=== extract feature from W ===")
    k = W.shape[1]
    m = W.shape[0]
    s_list = []
    for i in range(m):
        rowsum = np.sum(
            W[
                i,
            ]
        )
        p_iq_x_list = []
        for q in range(k):
            p_iq = W[i, q] / rowsum
            if p_iq != 0:
                tmp = p_iq * math.log(p_iq, 2)
            else:
                tmp = 0
            p_iq_x_list.append(tmp)
        s = 1 + 1 / math.log(k, 2) * np.sum(p_iq_x_list)
        s_list.append(s)
    return s_list


def predict_H(H):
    """extract feature from H"""
    print("=== extract feature from H ===")
    colmax = np.amax(H, axis=0)
    colsum = np.sum(H, axis=0)
    p = colmax / colsum
    idx = H.argmax(axis=0)
    out = [idx, p]
    return out


def cal_connectivity(H, idx):
    """calculate connectivity matrix"""
    print("=== calculate connectivity matrix ===")
    connectivity_mat = np.zeros((H.shape[1], H.shape[1]))
    classN = H.shape[0]
    for i in range(classN):
        xidx = list(np.concatenate(np.where(idx == i)))
        iterables = [xidx, xidx]
        for t in itertools.product(*iterables):
            connectivity_mat[t[0], t[1]] = 1
    return connectivity_mat


def cal_silhouette(C):
    silhouette = C
    return silhouette


def cal_cophenetic(C):
    """calculate cophenetic correlation coefficient"""
    print("=== calculate cophenetic correlation coefficient ===")
    X = C  # Original data (1000 observations)
    """Z = linkage(X)"""
    Z = fc.linkage_vector(X)  # Clustering
    orign_dists = fc.pdist(X)  # Matrix of original distances between observations
    cophe_dists = cophenet(Z)  # Matrix of cophenetic distances between observations
    corr_coef = np.corrcoef(orign_dists, cophe_dists)[0, 1]
    return corr_coef


def cal_dispersion(C):
    """calculate dispersion coefficient"""
    print("=== calculate dispersion coefficient ===")
    n = C.shape[1]
    corr_disp = np.sum(4 * np.square(np.concatenate(C - 1 / 2))) / (np.square(n))
    return corr_disp

def run_nmf(V, rank, n, prefix):
    """
    Run standard NMF on data set. n runs of Standard NMF are performed and obtained consensus matrix
    averages all n connectivity matrices.

    :param V: Target matrix with gene expression data.
    :type V: `sparse.matrix`
    :param rank: Factorization rank.
    :type rank: `int`
    :param n: # of runs
    """
    if n == 1:
        print("=== run NMF at rank %d  ===" % rank)
        model = NMF(n_components=rank, init="nndsvd", random_state=0, verbose=True)
        W = model.fit_transform(V)
        H = model.components_
        print(
            "1/%d : reconstruction err: %s (%3d/200 iterations)"
            % (n, model.reconstruction_err_, model.n_iter_)
        )
        """
        o_sparseH = cal_sparseness(H)
        o_sparseW = cal_sparseness(W)
        o_rss_mse = cal_rss_mse(W,H,V)
        o_rss = o_rss_mse[0]
        o_mse = o_rss_mse[1]
        o_evar = cal_evar(o_rss, V)
        """
        o_fsW = cal_featureScore_kim(W)
        o_predH = predict_H(H)
        """
        out = [rank, n, o_sparseH, o_sparseW, o_rss, o_mse, o_evar]
        np.savetxt('.'.join([prefix, "sta.txt"]), out)
        """
        np.savetxt(".".join([prefix, "featureScore_W.txt"]), o_fsW)
        np.savetxt(".".join([prefix, "predict_H.txt"]), o_predH)
        saveH(prefix, H)
        saveW(prefix, W)
    else:
        print("=== run NMF at rank %d with %d runs ===" % (rank, n))
        out_list = []
        consensus = np.zeros((V.shape[1], V.shape[1]))
        for i in range(n):
            model = NMF(
                n_components=rank,
                init="nndsvd",
                random_state=i,
                verbose=True,
                max_iter=500,
            )
            W = model.fit_transform(V)
            H = model.components_
            print(
                "%2d/%d : reconstruction err: %s (%3d/500 iterations)"
                % (i + 1, n, model.reconstruction_err_, model.n_iter_)
            )
            consensus += cal_connectivity(H, predict_H(H)[0])
            """
            o_sparseH = cal_sparseness(H)
            o_sparseW = cal_sparseness(W)
            o_rss_mse = cal_rss_mse(W,H,V)
            o_rss = o_rss_mse[0]
            o_mse = o_rss_mse[1]
            o_evar = cal_evar(o_rss, V)
            out = [i+1, rank, n, o_sparseH, o_sparseW, o_rss, o_mse, o_evar]
            out_list.append(out)
            """
        consensus /= n
        p_consensus = reorder(consensus)
        plotC(prefix, p_consensus, rank)
        saveC(prefix, p_consensus)
        """
        o_cophcor = cal_cophenetic(consensus)
        o_disp = cal_dispersion(consensus)
        np.savetxt('.'.join([prefix, "sta.mx"]), out_list, delimiter="\t")
        out2 = list(np.mean(np.squeeze(out_list)[:,1:], axis=0))
        out2.append(o_cophcor)
        out2.append(o_disp)
        np.savetxt('.'.join([prefix, "sta.txt"]), out2)
        """
        print("perform NMF in nndsvd model")
        model = NMF(n_components=rank, init="nndsvd", random_state=0, verbose=True)
        W = model.fit_transform(V)
        H = model.components_
        saveH(prefix, H)
        saveW(prefix, W)
        o_fsW = cal_featureScore_kim(W)
        o_predH = predict_H(H)
        np.savetxt(".".join([prefix, "featureScore_W.txt"]), o_fsW)
        np.savetxt(
            ".".join([prefix, "predict_H.txt"]), np.squeeze(o_predH).T, delimiter="\t"
        )
def run_nmf_2(V, rank, n, prefix):
    """
    Run standard NMF on data set. n runs of Standard NMF are performed and obtained consensus matrix
    averages all n connectivity matrices.

    :param V: Target matrix with gene expression data.
    :type V: `sparse.matrix`
    :param rank: Factorization rank.
    :type rank: `int`
    :param n: # of runs, which determines the random_state of nmf
    """
    print("perform NMF 2 in nndsvd model")
    model = NMF(
        n_components=rank,
        init="nndsvd",
        ## random_state based on n
        random_state=n,
        verbose=True,
        max_iter = 500
    )
    W = model.fit_transform(V)
    H = model.components_
    saveH(prefix, H)
    saveW(prefix, W)
    # o_fsW = cal_featureScore_kim(W)
    # o_predH = predict_H(H)
    # np.savetxt(".".join([prefix, "featureScore_W.txt"]), o_fsW)
    # np.savetxt(
    #     ".".join([prefix, "predict_H.txt"]), np.squeeze(o_predH).T, delimiter="\t"
    # )

def post_nmf_2(input_prefix,
               n_mods,
               out_prefix,
               start_random_state = 0,
               end_random_state = 9,
               default_random_state = 0):
    # load all the generated H
    H_candidates = [f"{input_prefix}.r{n_mods}.n{i}.H.mx"
             for i in range(start_random_state, end_random_state)]
    H_exists = [j for j in H_candidates if exists(j)]
    print("{len(H_exists)} H.mx files are found for r{n_mods}.")
    used_random_states = [basename(i).split(".")[3] for i in H_exists]
    print(f"Found random states: {','.join(used_random_states)}")
    # load default H and W
    H = np.loadtxt(
        f"{input_prefix}.r{n_mods}.n{default_random_state}.H.mx")
    W = np.loadtxt(
        f"{input_prefix}.r{n_mods}.n{default_random_state}.W.mx")
    # get consensus matrix (depend on Hs)
    n_cluster = H.shape[1]
    consensus = np.zeros(shape = (n_cluster, n_cluster))
    for ih in H_exists:
        iH = np.loadtxt(ih)
        consensus += cal_connectivity(iH, predict_H(iH)[0])
    consensus /= len(H_exists)
    # plot and save results
    p_consensus = reorder(consensus)
    outer = f"{out_prefix}.r{n_mods}"
    plotC(outer, p_consensus, n_mods)
    saveC(outer, p_consensus)
    o_fsW = cal_featureScore_kim(W)
    o_predH = predict_H(H)
    np.savetxt(".".join([outer, "featureScore_W.txt"]), o_fsW)
    np.savetxt(".".join([outer, "predict_H.txt"]),
               np.squeeze(o_predH).T, delimiter = "\t")
    


def plotC(prefix, C, rank):
    """
    Plot reordered consensus matrix.

    :param C: Reordered consensus matrix.
    :type C: numpy.ndarray`
    :param rank: Factorization rank.
    :type rank: `int`
    """
    fig = plt.figure(figsize=(5, 5), dpi=100)
    imshow(C)
    set_cmap("RdBu_r")
    fileN = [prefix, "C", "png"]
    fileN = ".".join(fileN)
    fig.savefig(fileN)


def reorder(C):
    """
    Reorder consensus matrix.

    :param C: Consensus matrix.
    :type C: `numpy.ndarray`
    """
    Y = 1 - C
    Z = linkage(squareform(Y), method="average")
    ivl = leaves_list(Z)
    ivl = ivl[::-1]
    return C[:, ivl][ivl, :]


if __name__ == "__main__":
    """Run NMF and save H & W."""
    run()
