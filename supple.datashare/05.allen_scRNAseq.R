library(tidyverse)
library(Seurat)

# meta
projd <- "/projects/ps-renlab2/szu/projects/CEMBA2"
allend <- file.path(projd, "19.snap2_integration", "out",
  "transferLabel_tscc")

# check Allen dataset
# this is all platform(10xv3, 10xv3, snRNA); male only;
# downsampled on subclass level (latest subclass)
r <- readRDS(file.path(allend, "allen_ds1000_seurat.rds"))
