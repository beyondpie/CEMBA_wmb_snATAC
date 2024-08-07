library(tidyverse)
library(Seurat)

# meta
projd <- "/projects/ps-renlab2/szu/projects/CEMBA2"
allend <- file.path(projd, "19.snap2_integration", "out",
  "transferLabel_tscc")

# check Allen dataset
r <- readRDS(file.path(allend, "allen_ds1000_seurat.rds"))
