library(tidyverse)
library(rlang)
library(R6)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

library(data.table)
library(GenomicRanges)
library(S4Vectors)

rdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env, Sa2PeakCalling)
import::from(.from = "peak.R", .directory = rdir,
  mapL4pc2L4)
work_dir <- file.path(here::here(), "18.snap2_peakcalling")

# * load union peaks
unionPeakBed <- Sa2PeakCalling$readUnionPeakBedFile()
clsUnionPeakBeds <- Sa2PeakCalling$loadpL4UnionPeakBeds()
# check
bedExists <- vapply(clsUnionPeakBeds,
  \(x) {
    all(rownames(x) %in% rownames(unionPeakBed))
  }, TRUE)

saveRDS(clsUnionPEakBeds,
  file = file.path(work_dir, "out/tscc/pL4UnionPeakBeds.rds"))

