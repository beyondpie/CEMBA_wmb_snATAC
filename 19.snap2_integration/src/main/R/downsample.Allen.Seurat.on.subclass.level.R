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
if(!dir.exists(outdir)) {
  dir.create(outdir)
}
allenAnnotMeta <- Sa2Integration$loadAllenAnnot()
allenAnnotMeta$subclass_label_v3 <- vapply(allenAnnotMeta$subclass_label,
  Sa2PeakCalling$renameAllenSubclass, "rename")

neuronAllenSeu <- readRDS(
  file.path(outdir, "neuron_allen_noraw_seurat.rds"))
neuronAllenSeu$subclass <- allenAnnotMeta[
  as.character(neuronAllenSeu$cl), "subclass_label_v3"]
neuronAllenSeu$barcode <- colnames(neuronAllenSeu)

nnAllenSeu <- readRDS(
  file.path(outdir, "nn_allen_noraw_seurat.rds"))
nnAllenSeu$subclass <- allenAnnotMeta[
  as.character(nnAllenSeu$cl), "subclass_label_v3"]
nnAllenSeu$barcode <- colnames(nnAllenSeu)

# * downsample by subclass
s <- merge(nnAllenSeu, neuronAllenSeu)
nds <- 1000
barcodes <- s@meta.data |>
  group_by(subclass) |> slice_sample(n = nds) |>
  x => x$barcode
s <- s[, barcodes]
mat1 <- as(s$RNA[["data.1"]], Class = "dgCMatrix")
mat2 <- as(s$RNA[["data.2"]], Class = "dgCMatrix")
mat <- SeuratObject::RowMergeSparseMatrices(
  mat1 = mat1, mat2 = mat2)
meta <- s@meta.data
r <- Seurat::CreateSeuratObject(counts = mat, assay = "RNA",
  meta.data = meta)
r <- Seurat::SetAssayData(object = r, slot = "data", new.data = mat)
saveRDS(r, file.path(outdir, "allen_ds1000_seurat.rds"))

## nds <- 5000
## neuron.meta.ds <- neuronAllenSeu@meta.data |>
##   group_by(subclass) |> slice_sample(n = nds)
## neuronAllenSeu.ds <- neuronAllenSeu[, neuron.meta.ds$barcode]
## saveRDS(neuronAllenSeu.ds, file.path(outdir,
##   "neuron_allen_ds5000_seurat.rds"))

## nn.meta.ds <- nnAllenSeu@meta.data |>
##   group_by(subclass) |> slice_sample(n = nds)
## nnAllenSeu.ds <- nnAllenSeu[ , nn.meta.ds$barcode]
## saveRDS(nnAllenSeu.ds, file.path(outdir,
##   "nn_allen_ds5000_seurat.rds"))

