import os
import sys
import traceback
from pathlib import Path
from typing import Dict

import pandas as pd
import snapatac2 as sa2
from snapatac2 import AnnData
code_root_dir = "/projects/ps-renlab/szu/projects/CEMBA2"
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils #pyright: ignore # noqa: F401, E402

# * log
logger = utils.set_file_logger( #pyright: ignore
    fnm = snakemake.log[0], #pyright: ignore # noqa: F821
    name = "sa2.embed"
)
# logger = utils.set_file_logger( #pyright: ignore
#     fnm = "test_qc_dlt.log", #pyright: ignore # noqa: F821
#     name = "sa2.embed"
# )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

def get_sample2bamfile(
        f: str = f"{code_root_dir}/meta/sample2rawbam.csv"):
    s2b = pd.read_csv(f, header = None, names = ['sample', 'bamfile', 'tag'])
    r = {row['sample']: row['bamfile'] for _, row in s2b.iterrows()}
    return r

def get_chromsize(f:str = f"{code_root_dir}/meta/mm10.chrom.sizes.lite"):
    chr2l = pd.read_csv(f, header = None, names = ['chr', 'len'],
                        delim_whitespace=True)
    r = {row['chr']: row['len'] for _, row in chr2l.iterrows()}
    return r
       
def QC_from_bamfile(chrom_sizes: Dict[str, int],
                    bam_file: str,
                    frag_file:str,
                    frag_stat_file: str,
                    snap_file: str,
                    jump_frag_file_if_exist: bool = False) -> AnnData:
    if not (os.path.exists(frag_file) and jump_frag_file_if_exist):
        stat = sa2.pp.make_fragment_file(
            bam_file = Path(bam_file),
            output_file = Path(frag_file),
            is_paired = True,
            barcode_regex = "^(\w+):.+", #pyright: ignore
            shift_left = 4,
            shift_right = -5,
            min_mapq = 30
        )
        with open(frag_stat_file, 'w') as f:
            f.writelines(str(stat))
    snap = sa2.pp.import_data(
        fragment_file = Path(frag_file),
        genome = sa2.genome.mm10, # type: ignore
        file = Path(snap_file),
        min_num_fragments = 200,
        min_tsse = 1,
        shift_left = 0,
        shift_right = 0,
        chrom_size = chrom_sizes,
        sorted_by_barcode = False)
    sa2.pp.filter_cells(
        data = snap,
        min_counts = 1000,
        min_tsse = 10,
        inplace = True)
    return snap

def DoubletRemoval_with_bmat(snap: AnnData, blacklist_file,
                             n_features = 99999999) -> None:
    """
    Use snapatac2 to
    1. transform original bam files to fragment files.
    2. for each fragment file, we perform qc
    - #fragment/cell >= 1000
    - tsse/cell >= 10
    3. then use snapatac2 to do doublet removal
    - instead of using gmat as what we did before,
        we use bmat after feature selection.
    - if we use all feature, it takes at least 12 hours.
    """
    sa2.pp.add_tile_matrix(
        adata = snap,
        bin_size = 500,
        inplace = True,
        exclude_chroms = ['chrM', 'M', 'chrUn*']
    )
    sa2.pp.select_features(snap,
                           n_features = n_features,
                           filter_lower_quantile = 0.005,
                           filter_upper_quantile = 0.005,
                           blacklist = blacklist_file,
                           max_iter = 1,
                           inplace = True)
    sa2.pp.scrublet(adata = snap,
                    features = 'selected',
                    expected_doublet_rate = 0.08,
                    random_state = 0,
                    inplace = True)
    return None
    
# if __name__ == '__main__':
#     # sample = sys.argv[0]
#     # out_dir = "sa2_QC_DoubletRemoval"
#     # test
#     # needs about 1 hour
#     sample = "CEMBA171206_3C"
#     out_dir = "test_qc_dlt"
#     os.makedirs(out_dir, exist_ok = True)
#     sample2bamfile = get_sample2bamfile()
#     chrom_sizes = get_chromsize()
#     blacklist_file = f"{code_root_dir}/meta/mm10.blacklist.bed"
#     snap_after_qc = QC_from_bamfile(
#         chrom_sizes = chrom_sizes, # pyright: ignore
#         bam_file = str(sample2bamfile[sample]),
#         frag_file = f"{out_dir}/{sample}_frag.h5ad",
#         frag_stat_file = f"{out_dir}/{sample}_frag_stat.txt",
#         snap_file = f"{out_dir}/{sample}_sa2_qc_dlt.h5ad"
#     )
#     snap_after_dlt = DoubletRemoval_with_bmat(
#         snap = snap_after_qc,
#         blacklist_file = blacklist_file
#     )
#     snap_after_dlt.close()

    
sample2bam_file = snakemake.input['sample2bam_file'] #pyright: ignore # noqa: F821
chromsizes_file = snakemake.input['chromsizes_file'] #pyright: ignore # noqa: F821
blacklist_file = snakemake.input['blacklist_file'] #pyright: ignore # noqa: F821

frag_file = snakemake.output['frag_file'][0] #pyright: ignore # noqa: F821
frag_sum_file = snakemake.output['frag_sum_file'][0] #pyright: ignore # noqa: F821
qc_dlt_file = snakemake.output['qc_dlt_file'][0] #pyright: ignore # noqa: F821

sample = snakemake.params['sample'] #pyright: ignore # noqa: F821

# sample = "CEMBA171206_3C"
# code_dir = "/projects/ps-renlab/szu/projects/CEMBA2"
# out_dir = "/oasis/tscc/scratch/szu/projects/CEMBA2/17.snapatac2/sa2_qc_dlt"
# sample2bam_file = f"{code_dir}/meta/sample2rawbam.csv"
# chromsizes_file = f"{code_dir}/meta/mm10.chrom.sizes.lite"
# blacklist_file = f"{code_dir}/meta/mm10.blacklist.bed"
# frag_dir = f"{out_dir}/frag"
# frag_sum_dir = f"{out_dir}/frag_sum"
# qc_dlt_dir = f"{out_dir}/qc_dlt"
# os.makedirs(frag_dir, exist_ok = True)
# os.makedirs(frag_sum_dir, exist_ok = True)
# os.makedirs(qc_dlt_dir, exist_ok = True)
# frag_file = f"{frag_dir}/{sample}_frag.h5ad"
# frag_sum_file = f"{frag_sum_dir}/{sample}_frag_sum.txt"
# qc_dlt_file = f"{qc_dlt_dir}/{sample}_qc_dlt.h5ad"

sample2bamfile = get_sample2bamfile(sample2bam_file)
chrom_sizes = get_chromsize(chromsizes_file)
logger.info(f"QC for {sample} to {frag_file}.")
s_qc = QC_from_bamfile(
    chrom_sizes = chrom_sizes, #pyright: ignore
    bam_file = str(sample2bamfile[sample]),
    frag_file = frag_file,
    frag_stat_file = frag_sum_file,
    snap_file = qc_dlt_file
)
logger.info(f"Doublet removal for {sample} to {qc_dlt_file}.")
DoubletRemoval_with_bmat(
    snap = s_qc,
    blacklist_file = blacklist_file
)
s_qc.close()
logger.info(f"QC and DoubletRemoval are done for {sample}.")
