library(BPCells)
library(Seurat)
library(SeuratObject)
options(Seurat.object.assay.version = "v5")
# Reference:
# https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette

library(future)
plan(multicore, workers = 2)
options(future.globals.maxSize = 20000 * 1024^2)
options(future.rng.onMisuse="ignore")


library(reticulate)
reticulate::use_python("/home/szu/mambaforge/envs/sa2/bin/python")
ad <- reticulate::import("anndata", convert = FALSE)

library(ggplot2)
library(purrr)
library(furrr)
library(dplyr)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)


## library(SeuratObject)
## library(SeuratDisk)

proj_root <- here::here()
rsc_dir <- "../resource/"
rsc2_dir <- file.path(proj_root, "19.snap2_integration",
  "src/main/resource")


## ann_fnm <- file.path(rsc_dir, "nn_male_10xv3_ann_noraw.h5ad")
## xlognorm <- BPCells::open_matrix_anndata_hdf5(path = ann_fnm, group = "X")
## BPCells::write_matrix_dir(
##   mat = xlognorm,
##   dir = file.path(rsc_dir, "nn_10xv3_ann")
## )

## norm_data <- BPCells::open_matrix_dir(file.path(rsc_dir, "nn_10xv3_ann"))

## obj1 <- Seurat::CreateSeuratObject(counts = norm_data)
## obj1 <- SetAssayData(object = obj1,
##   slot = 'data', new.data =  norm_data)

## ann <- ad$read_h5ad(ann_fnm, backed = 'r')
## obs <- py_to_r(ann$obs)
## # TODO: as.characters to convert factor to string

## # now let's test integration ability

## obj2 <- Seurat::CreateSeuratObject(
##   counts = BPCells::open_matrix_dir(
##     file.path(rsc_dir, "nn_10xv3_ann_bk")))
 
## obj2 <- SetAssayData(object = obj2,
##   slot = 'data', new.data =  norm_data)

## obj1 <- FindVariableFeatures(obj1, selection.method = "vst",
##   nfeatures = 2000)
## features <- VariableFeatures(obj1)
## obj1$group <- "g1"
## obj2$group <- "g2"

## obj1 <- ScaleData(obj1)
## obj2 <- ScaleData(obj2)

## # which takes about 30 mins
## # 23.7G for 509,173 cells x  2,000 features
## obj1 <- RunPCA(obj1, npcs = 50, features = features)
## obj2@reductions$pca <- obj1@reductions$pca

## # - computing within dataset neighborhoods
## # two extra Rs (so proportion to object.list)
## # each cost 16G
## # - when finding all pairwise anchors
## # only one extra R.
## # costs 11.6G
## # rpca is fast
## anchors <- FindIntegrationAnchors(
##   object.list = list(obj1, obj2),
##   normalization.method = NULL,
##   anchor.features = features,
##   scale = FALSE,
##   reduction = "rpca",
##   l2.norm = TRUE,
##   dims = 1:30,
##   verbose = TRUE
## )

## anchors.cca <- IntegrateLayers(
##   object = obj1,
##   method = CCAIntegration,
##   orig.reduction = "pca",
##   new.reduction = "integratied.cca",
##   verbose = TRUE,
##   features = features
## )

## ## CCA now works under IntegrateLayers
## # but it needs pca, and I guess it
## # perform CCA on PCA, like Harmony, RPCA.
## # but takes lots of time, so stop it
## merged_object <- merge(
##   x = obj1, y = obj2, add.cell.ids = c("Obj1", "Ob2"),
##   merge.data = TRUE,
##   merge.dr = TRUE)
## anchors.cca <- IntegrateLayers(
##   object = merged_object,
##   method = CCAIntegration,
##   orig.reduction = "pca",
##   new.reduction = "integratied.cca",
##   verbose = TRUE,
##   features = features
## )

## # Error: S4, but double found
## anchors.cca <- FindIntegrationAnchors(
##   object.list = list(obj1, obj2),
##   normalization.method = "LogNormalize",
##   anchor.features = features,
##   scale = TRUE,
##   reduction = "cca",
##   l2.norm = TRUE,
##   dims = 1:30,
##   verbose = TRUE
## )

## # * now test neurons
## BPCells::write_matrix_dir(
##   mat = BPCells::open_matrix_anndata_hdf5(
##     path = file.path(rsc2_dir,
##       "norawdata_allen",
##       "neuron_male_10xv3_ann_noraw.h5ad"),
##     group = "X"),
##   dir = file.path(rsc_dir, "neuron_10xv3")
## )

## BPCells::write_matrix_dir(
##   mat = BPCells::open_matrix_anndata_hdf5(
##     path = file.path(rsc2_dir,
##       "norawdata_atac",
##       "neuron_gmat_atac_ann.h5ad"),
##     group = "X"
##   ),
##   dir = file.path(rsc_dir, "neuron_atac")
## )

allen <- CreateSeuratObject(
  counts = BPCells::open_matrix_dir(file.path(rsc_dir, "neuron_10xv3")),
  assay = "RNA")

allen <- SetAssayData(object = allen, slot = "data",
  new.data = BPCells::open_matrix_dir(file.path(rsc_dir, "neuron_10xv3"))
)
  
atac <- CreateSeuratObject(
  counts = BPCells::open_matrix_dir(file.path(rsc_dir, "neuron_atac")),
  assay = "ATAC") |>
  SetAssayData(object = _, slot = "data",
    new.data = BPCells::open_matrix_dir(file.path(rsc_dir, "neuron_atac")))

allen <- FindVariableFeatures(allen, selection.method = "vst",
  nfeatures = 2000)
features <- VariableFeatures(allen)

allen$group <- "allen"
atac$group <- "atac"

# 19G for each one
allen_atac <- furrr::future_map(list(allen, atac), \(x) {
  ScaleData(x, features = features)
})

# one worker cost 50G RAM (should be atac)
# master and the other work costs 20G RAM
# needs about 2 hours
allen_atac <- furrr::future_map(allen_atac, \(x) {
  RunPCA(x, npcs = 50, features = features)
})

saveRDS(allen_atac, "tmp.allen_atac.pca.rds")

allen_atac <- readRDS("tmp.allen_atac.pca.rds")
features <- VariableFeatures(allen_atac[[1]])

anchors.neuron <- FindIntegrationAnchors(
  object.list = allen_atac,
  normalization.method = NULL,
  anchor.features = features,
  scale = FALSE,
  reduction = "rpca",
  l2.norm = TRUE,
  dims = 1:30,
  verbose = TRUE
)
saveRDS(anchors.neuron, "tmp.anchors.neuron.rds")


## check FindIntegrationAnchors
merged_allen_atac <- merge(
  x = allen_atac[[1]],
  y = allen_atac[[2]],
  add.cell.ids = NULL,
  merge.data = TRUE,
  merge.dr = TRUE)

## and IntegrationAnalysis on rpca
anchors.rpca <- IntegrateLayers(
  object = merged_allen_atac,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integratied.rpca",
  verbose = TRUE,
  features = features,
  scale.layer = "scale.data",
  scale = FALSE,
  normalization.methpd = "LogNormalize",
  group.by = "group"
)


# 32G at least
anchors.cca <- IntegrateLayers(
  object = merged_allen_atac,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integratied.rpca",
  verbose = TRUE,
  features = features,
  scale.layer = "scale.data",
  scale = FALSE,
  normalization.methpd = "LogNormalize",
  group.by = "group"
)
