library(tidyverse)
library(Seurat)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 5e9)
library(SeuratObject)
library(Matrix)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(future)

projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cluSumBySa2, Sa2Integration, cembav2env, Sa2PeakCalling)
import::from(.from = "integration.R", .directory = rdir,
  convertAnn2Seurat5, get.downsample.fun, downsampleSeurat,
  toSeuratInMemory,
  isOnDiskMat.Seurat, calVarOfFea.Seurat, setVariableFeatures)

# * configs
# tscc or encoder
system <- "tscc"
rscdir <- file.path(here::here(), "19.snap2_integration",
  "src/main/resource")
outdir <- file.path(here::here(), "19.snap2_integration",
  paste0("out/transferLabel_", system))

# load atacmeta
atacMeta <- readRDS(cemabv)
