from typing import List
import os
import snapatac2 as sa2
import pyprojroot

proj_root = pyprojroot.here()
raw_bam_dir: str = os.path.join(proj_root, "data",
                                "raw_bam_test")
dedup_bam_dir: str = os.path.join(proj_root, "data",
                                  "r")
samples: List[str] = [
    "CEMBA181023_6B", "CEMBA201210_10D"]
outdir = os.path.join(proj_root, "00.data.preprocess",
                      "out/test")

for s in samples:
    stat = sa2.pp.make_fragment_file(
        bam_file = os.path.join(raw_bam_dir, f"{s}.bam"),
        output_file = os.path.join(outdir,
                                   f"{s}.frag.rawbam.tsv"),
        is_paired = True,
        barcode_regex = "^(\w+):.+",
        shift_left = 4,
        shift_right = -4,
        min_mapq = 30
    )
    print(stat)

for s in samples:
    stat = sa2.pp.make_fragment_file(
        bam_file = os.path.join(proj_root, "data",
                                "filtered_dedup_bam",
                                f"{s}.filtered_dedup.bam"),
        output_file = os.path.join(outdir,
                                   f"{s}.frag.filtered_debup_bam.tsv"),
        is_paired = True,
        barcode_regex = "^(\w+):.+",
        shift_left = 4,
        shift_right = -4,
        min_mapq = 30
    )

