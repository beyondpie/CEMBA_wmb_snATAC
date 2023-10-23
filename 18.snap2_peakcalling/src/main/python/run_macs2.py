import os
import sys
import logging
from dataclasses import dataclass, field
from typing import Any, List, Dict
from subprocess import Popen, PIPE, STDOUT, check_output

import pyprojroot
proj_dir = str(pyprojroot.here())
py_dir = f"{proj_dir}/package/python"
sys.path.insert(0, py_dir)
from mylog import StreamToLogger, set_file_logger

# * snakemake
logfnm: str = snakemake.log[0]
bedfile_dirs: Dict[str, str] = snakemake.params["bedfile_dirs"]
outdir: str = snakemake.params["outdir"]
prefix: str = snakemake.params["prefix"]
cl: str = snakemake.wildcards["cluster"]
nmin: int = snakemake.params["nmin"]
macs2: str = snakemake.params["macs2"]
debug: int = snakemake.params["debug"]

# * debug
if debug > 0:
    print("Under debug mode ...")
    system: str = "imac"
    out_dir: str = os.path.join(proj_dir, "18.snap2_peakcalling/out")
    logfnm: str = os.path.join(out_dir, "log", "test_macs2.log")
    outdir: str = os.path.join(out_dir, system, "macs2")
    bedfile_dirs: Dict[str, str] = {
        "all": os.path.join(outdir, "all_bed"),
        "biorep": os.path.join(outdir, "biorep_bed"),
        "pseudorep": os.path.join(outdir, "pseudorep_bed")
    }
    macs_out_dir: str = os.path.join(outdir, "macs2")
    prefix: str = "neuron."
    cl: str = "2-1-1-16"
    size: List[int] = [714, 370, 344]
    nmin: int = 200
    macs2:str = "/home/szu/mambaforge/envs/sa2/bin/macs2"
else:
    print("Under normal mode ...")


# * set logger
logger = set_file_logger(logfnm, name = "run_macs2")
## works in Linux, but have some troubles in mac
sys.stdout = StreamToLogger(logger = logger, level = logging.INFO)
sys.stderr = StreamToLogger(logger = logger, level = logging.ERROR)

# * class and functions
def get_bedfile(bedir:str, prefix:str,
                cluster:str,
                suffix:str = ".bed.gz",
                check_exist: bool = True) -> str | None:
    """
    bedfile name is based on snapatac2.ex.export_bed logic.
    """
    r = str(os.path.join(bedir, f"{prefix}{cluster}{suffix}"))
    if not os.path.exists(r):
        if not check_exist:
            print(f"{r} does not exist, return None.")
            return None
        else:
            return r
    return r

def get_macs2_exp(macs2: str, bedfnm, outdir: str, name: str,
                  shift: int = -75, extsize: int = 150,
                  use_str: bool = True) -> str | List[str]:
    """
    FIXME: use_str is not recommended in python subprocess.
    But it does not work when I set it as False.
    - It cannot recognize ~ as home
    - It cannot recognize the parameters for macs2
    While set it as True and use shell as True then, it works.
    """
    r = [macs2, "callpeak", f"-t {bedfnm}",
         "-f BED",
         f"-n {name}", "-g mm",
         "-q 0.01", "--nomodel", f"--shift {shift}",
         f"--extsize {extsize}", "-B",
         "--SPMR", "--keep-dup all", "--call-summits",
         f"--outdir {outdir}"]
    if use_str:
        return ' '.join(r)
    # Popen prefer this
    else:
        return r
    
def log_subprocess_output(pipe, logger, prefix:str ):
    """
    ref: https://stackoverflow.com/questions/21953835/run-subprocess-and-print-output-to-logging
    """
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        logger.info(f"{prefix}: {line}")
        
@dataclass
class Cluster4Macs2:
    cluster:str
    all_bedfile: str|None= None
    early_bedfile: str|None = None
    later_bedfile: str|None = None
    earlypseudo_bedfile: str|None = None
    laterpseudo_bedfile: str|None = None

    def __init__(self, cluster: str):
        self.cluster = cluster

    def set_bedfile(self,
                    prefix: str,
                    all_bedir: str,
                    biorep_bedir: str,
                    pseudorep_bedir: str,
                    suffix: str = ".bed.gz") -> None:
        self.all_bedfile = get_bedfile(
            all_bedir, prefix, self.cluster, suffix, True)
        self.early_bedfile = get_bedfile(biorep_bedir, prefix,
                                         f"{self.cluster}.early", suffix, False)
        self.later_bedfile = get_bedfile(biorep_bedir, prefix,
                                         f"{self.cluster}.later", suffix, False)
        self.earlypseudo_bedfile = get_bedfile(pseudorep_bedir, prefix,
                                               f"{self.cluster}.early.pseudo",
                                               suffix, False)
        self.laterpseudo_bedfile = get_bedfile(pseudorep_bedir, prefix,
                                               f"{self.cluster}.later.pseudo",
                                               suffix, False)
        
# * main
logger.info(f"Initilize cluster of {cl} for macs2 info.")
cl4macs2 = Cluster4Macs2(cluster = cl)
logger.info(f"set the bedfiles for cluster {cl}.")
cl4macs2.set_bedfile(prefix = prefix,
                     all_bedir = bedfile_dirs["all"],
                     biorep_bedir=bedfile_dirs["biorep"],
                     pseudorep_bedir=bedfile_dirs["pseudorep"])
use_str_in_shell: bool = True
process_macs2_pool = []

for i in [cl4macs2.all_bedfile,
          cl4macs2.early_bedfile, cl4macs2.later_bedfile,
          cl4macs2.earlypseudo_bedfile, cl4macs2.laterpseudo_bedfile]:
    if i is not None:
        nm = os.path.basename(i).rstrip(".bed.gz")
        macs2_rep = get_macs2_exp(macs2 = macs2,
                                       bedfnm = i,
                                       outdir = outdir,
                                       name = nm,
                                       use_str = use_str_in_shell)
        logger.info(f"start macs2 for : {nm} with command: {macs2_rep}")
        process_macs2_pool.append(Popen(macs2_rep,
                                        shell = use_str_in_shell))

exit_signals = [p.wait() for p in process_macs2_pool]
if sum(exit_signals) != 0:
    raise RuntimeError(
        f"macs2 returns {','.join(exit_signals)} instead of all zeros.")
else:
    logger.info("all macs2 process done.")
