import sys
import traceback
from pathlib import Path

import snapatac2 as sa2

import pyprojroot
code_root_dir = str(pyprojroot.here())
pack_dir = f"{code_root_dir}/package/python"
sys.path.insert(0, pack_dir)
import utils #pyright: ignore # noqa: F401, E402

# * log
logger = utils.set_file_logger( #pyright: ignore
    fnm = snakemake.log[0], #pyright: ignore # noqa: F821
    name = "sa2.get.sample.gmat"
)
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(
                             exc_type, exc_value, exc_traceback)
                         ]))
# Install exception handler
sys.excepthook = handle_exception

snap_file = snakemake.input["snap_file"][0]
genome = snakemake.params['genome']
out_file = snakemake.output["gmat_file"][0]

logger.info(f"Load snap file {snap_file}")
snap = sa2.read(snap_file, backed = 'r')
# NOTE:
# currently, only sa2 default genome support (i.e., mm10 in mouse)
logger.info(f"genome: {genome}.")
logger.info(f"write gmat to file: {out_file}.")
sa2.pp.make_gene_matrix(
    adata = snap,
    gene_anno = sa2.genome.mm10,
    file = Path(out_file),
    use_x = False,
    id_type = "gene"
)

logger.info("Done.")
