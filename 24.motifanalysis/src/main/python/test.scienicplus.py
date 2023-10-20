import pycistarget
import pyranges as pr
import pandas as pd
from scenicplus.wrappers.run_pycistarget import run_pycistarget
import glob
import os
from pyscistarget.utils import region_names_to_coordinates

import pyprojroot

proj_dir = str(pyprojroot.here())

nvolp_peak_fnm = os.path.join(proj_dir, "18.snap2_peakcalling",
                              "out/scfilter", "cembav2.nonOvlpDHS.bed")

peaks = pd.read_csv(nvolp_peak_fnm, sep = "\t", header = None)
peaks = peaks.rename(columns = {0: "Chromosome", 1: "Start", 2: "End", 3: "Name"})

region_sets = {'DARs':
               {'novlp': pr.PyRanges(peaks[["Chromosome", "Start", "End"]])}
               }
db_path = os.path.join(proj_dir, "24.scenicplus/data/scenicplus_database")
rankings_db = os.path.join(db_path, "mm10_screen_v10_clust.regions_vs_motifs.rankings.feather")
scores_db = os.path.join(db_path, "mm10_screen_v10_clust.regions_vs_motifs.scores.feather")
motif_annotation = os.path.join(db_path, "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")

out_dir = os.path.join(proj_dir, "24.sceinicplust/out/test_scenicplus")

motif_res = run_pycistarget(
    region_sets = region_sets,
    species = "mus_musculus",
    save_path = out_dir,
    ctx_db_path = rankings_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    annotation_version = "v10nr_clust",
    n_cpu = 4,
    ignore_reinit_error = True,
)
