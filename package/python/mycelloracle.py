import pandas as pd
import numpy as np
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm
import celloracle as co
from celloracle import motif_analysis as ma
from genomepy import Genome
import scanpy as sc

def download_refgenome(ref = "mm10") -> None:
    genome_installation = ma.is_genome_installed(ref_genome = ref)
    if not genome_installation:
        import genomepy
        ## fasta: /Users/szu/.local/share/genomes/mm10/mm10.fa
        genomepy.install_genome(name = ref, provider = "UCSC")
    else:
        print(f"{ref} is installed")


def decompose_chrstr(peak_str):
    """
    Args:
        peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

    Returns:
        tuple: chromosome name, start position, end position
    """
    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)
    return chr_, start, end

def check_peak_format(peaks_df, ref_genome):
    """
    Check peak format.
     (1) Check chromosome name.
     (2) Check peak size (length) and remove sort DNA sequences (<5bp)

    """
    df = peaks_df.copy()
    n_peaks_before = df.shape[0]
    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed), index=peaks_df.index)
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(int)
    df_decomposed["end"] = df_decomposed["end"].astype(int)
    # Load genome data
    genome_data = Genome(ref_genome)
    all_chr_list = list(genome_data.keys())
    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])
    # Filter peaks with invalid chromosome name
    n_threshold = 5
    df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]
    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

    # Data counting
    n_invalid_length = len(lengths[lengths < n_threshold])
    n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
    n_peaks_after = df.shape[0]
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)
    return df

def load_tfidf(f) -> pd.DataFrame:
    import pyarrow.parquet as pq
    r = pq.read_table(f)
    return r.to_pandas()
