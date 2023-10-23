## This is for preprocess seurat v5 data
library(tidyverse)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratObject)
# Not sure if needed
library(SeuratData)
# seurat5 version wrappers
# needed for FastMNNIntegration
library(SeuratWrappers)
# [FIXED] Seurat5 use future inside, this is a must
options(future.globals.maxSize = 5e9)
# other packages
library(stringr)
library(purrr)
library(Matrix)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

library(reticulate)
options(stringsAsFactors = FALSE)
usePython <- "/home/szu/mambaforge/envs/sa2/bin/python"
reticulate::use_python("/home/szu/mambaforge/envs/sa2/bin/python")

projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cluSumBySa2, Sa2Integration, cembav2env)
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
geneList <- Sa2Integration$getMarkersList()
allenAnnotMeta <- Sa2Integration$loadAllenAnnot()

# ===================================================
# * transform allen's ann to Seurat
# add subclass, supertype, cluster_id to allenS5
toSeuratFromAllen <- function(cellfrom, tech,
                              meta = allenAnnotMeta) {
  annfnm <- file.path(rscdir, "norawdata_allen",
    str_glue("{cellfrom}_male_{tech}_ann_noraw.h5ad"))
  message("Read annfnm: ", annfnm)
  allenSeu <- convertAnn2Seurat5(
    annfnm = annfnm,
    modality = "rna",
    group = "X",
    outdir = file.path(outdir, str_glue("{cellfrom}_{tech}_noraw")),
    isLogNorm = TRUE,
    removeCounts = FALSE,
    overwrite = FALSE,
    addObsMeta = TRUE,
    usePython = "/home/szu/mambaforge/envs/sa2/bin/python"
  )
  allen_cl <- as.character(allenSeu$cl)
  allenSeu$class_id <- meta[allen_cl, "class_id"]
  allenSeu$subclass_id <- meta[allen_cl, "subclass_id"]
  allenSeu$supertype_id <- meta[allen_cl, "supertype_id"]
  outfnm <- file.path(outdir,
    str_glue("{cellfrom}_{tech}_noraw_seurat.rds"))
  message("save seurat to: ", outfnm)
  saveRDS(allenSeu, outfnm)
  return(allenSeu)
}

neuronAllenSeu <- toSeuratFromAllen(cellfrom = "neuron",
  tech = "allen")
nnAllenSeu <- toSeuratFromAllen(cellfrom = "nn",
  tech = "allen")
neuron10xv3Seu <- toSeuratFromAllen(cellfrom = "neuron",
  tech = "10xv3")
nn10xv3Seu <- toSeuratFromAllen(cellfrom = "nn",
  tech = "10xv3")

# ** transform atac's ann to Seurat
barcode2L4 <- cluSumBySa2$loadbarcode2L4()
toSeuratFromSa2ATAC <- function(cellfrom) {
  annfnm <- file.path(rscdir, "norawdata_atac",
    str_glue("{cellfrom}_gmat_atac_ann.h5ad"))
  atacSeu <- convertAnn2Seurat5(
    annfnm = annfnm,
    modality = "atac",
    group = "X",
    outdir = file.path(outdir, str_glue("{cellfrom}_atac_noraw")),
    isLogNorm = TRUE,
    removeCounts = FALSE,
    overwrite = FALSE,
    addObsMeta = FALSE,
    usePython = ""
  )
  atacSeu$L4 <- barcode2L4[colnames(atacSeu), "L4"]
  outfnm <- file.path(outdir,
    str_glue("{cellfrom}_atac_noraw_seurat.rds"))
  message("save to: ", outfnm)
  saveRDS(atacSeu, outfnm)
  return(atacSeu)
}
neuatacSeu <- toSeuratFromSa2ATAC(cellfrom = "neuron")
nnatacSeu <- toSeuratFromSa2ATAC(cellfrom = "nn")

# =========================================================
# * prepare allen data related with imn
# show columns
allenSubclass2Region <- Sa2Integration$getAllenSubclass2Region()
as.tibble(colnames(allenSubclass2Region)) |> print(10)
imnscs <- unique(with(allenAnnotMeta,
  subclass_id_label[grepl("IMN", subclass_id_label)]))
# [1] "038 DG-PIR Ex IMN"      "045 OB-STR-CTX Inh IMN"
n_imnscs_v3 <- unique(with(allenAnnotMeta,
  v3.size[subclass_label %in% imnscs]))
n_imnscs_v2 <- unique(with(allenAnnotMeta,
  v2.size[subclass_label %in% imnscs]))
# > 6000 cells
n_imnscs_rough <- (
  sum(n_imnscs_v3) + sum(n_imnscs_v2)) / 2

# get major regions of IMN
imn_region <- allenSubclass2Region[imnscs, "MajorRegionRateTop3"] |>
  lapply(X = _, Sa2Integration$extractAllenRegionWithScore) |>
  do.call(rbind, args = _) |>
  group_by(region) |>
  summarise(mean = sum(score))

imn_neuron_scs <- allenSubclass2Region$subclass_id_label |>
  vapply(X = _, \(r) {
    if (grepl(" NN$", r)) {
      return(FALSE)
    }
    rs <- allenSubclass2Region[r, "MajorRegionRateTop3"] |>
      Sa2Integration$extractAllenRegionWithScore()
    if (any(rs$region %in% imn_region$region )) {
      return(TRUE)
    }
    return(FALSE)
  }, FUN.VALUE = TRUE) |>
  x => allenSubclass2Region$subclass_id_label[x]

# about 1.1 million cells in total
imn_neuron_cls <- allenAnnotMeta$cl[allenAnnotMeta$subclass_id_label %in%
                                      imn_neuron_scs]
write.table(imn_neuron_cls, file = file.path(
  projdir, "19.snap2_integration", "src/main/resource",
  "AllenIMNTopMajorRegionRelatedAllenCls.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

# get L4 in these regions
imn_cemba_region <- allenSubclass2Region[imnscs, "MajorRegionRateTop3"] |>
  lapply(X = _, Sa2Integration$extractAllenRegionWithScore) |>
  do.call(rbind, args = _) |>
  group_by(region) |>
  summarise(mean = sum(score))

atacMeta <- readRDS(cembav2env$sa2metaFile)
a <- unique(atacMeta[, c("L4", "mainclass.rough", "nt_type_label.allen")])
rownames(a) <- a$L4
## write.table(a,
##   file = file.path(projdir, "19.snap2_integration", "src/main/resource",
##     "L4toMainClasstoAllenNtType.csv"),
##   sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

L4toRegion <- cluSumBySa2$loadL4MajorRegion()
imn_neuron_L4s <- L4toRegion$L4 |>
  vapply(X = _, \(r) {
    if (a[r, "mainclass.rough"] == "NN") {
      return(FALSE)
    }
    rs <- L4toRegion[r, "MajorRegionRateTop3"] |>
      Sa2Integration$extractAllenRegionWithScore()
    if (any(rs$region %in% imn_region$region)) {
      return(TRUE)
    }
    return(FALSE)
  }, FUN.VALUE = TRUE) |>
  x => L4toRegion$L4[x]
write.table(imn_neuron_L4s, file = file.path(
  projdir, "19.snap2_integration", "src/main/resource",
  "AllenIMNTopMajorRegionRelatedL4s.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

# ** now get the seurat data for them
imn_neuron_cls <- data.table::fread(
  file = file.path(rscdir, "AllenIMNTopMajorRegionRelatedAllenCls.txt"),
  header = FALSE, data.table = FALSE)$V1
neuronAllenSeu <- readRDS(
  file.path(outdir, "neuron_allen_noraw_seurat.rds"))
nnAllenSeu <- readRDS(
  file.path(outdir, "nn_allen_noraw_seurat.rds")
)

imnNeuronAllenSeu <- subset(neuronAllenSeu, subset = cl %in% imn_neuron_cls)
imnNNAllenSeu <- subset(nnAllenSeu, subset = cl %in% imn_neuron_cls)
imnNNeuronAllenSeu <- merge(imnNeuronAllenSeu, imnNNAllenSeu)
saveRDS(imnNNeuronAllenSeu,
  file = file.path(outdir, "imneuron_nneuronallen_noraw_seurat.rds"))

# Above on have issues to memory
# let's create on with downsampled directly
imnNNeuronAllenSeu <- readRDS(
  file.path(outdir, "imneuron_nneuronallen_noraw_seurat.rds"))
dsimnSeu <- downsampleSeurat(imnNNeuronAllenSeu, groupBy = "cl",
  minNum = 50, maxNum = 51)

mat1 <- as(dsimnSeu$RNA[["data.1"]], Class = "dgCMatrix")
mat2 <- as(dsimnSeu[["RNA"]][["data.2"]], Class = "dgCMatrix")
mat <- SeuratObject::RowMergeSparseMatrices(mat1 = mat1, mat2 = mat2)
meta <- dsimnSeu@meta.data
# all(rownames(meta) == colnames(mat))
r <- Seurat::CreateSeuratObject(counts = mat, assay = "RNA",
  meta.data = meta)
r <- Seurat::SetAssayData(object = r, slot = "data", new.data = mat)
saveRDS(r, file.path(outdir, "imneuron_dp50allen_noraw_seurat.rds"))

## saveRDS(imnNeuronAllenSeu,
##   file = file.path(outdir, "imneuron_allen_noraw_seurat.rds"))
## ** Allen's IMN
imn_neuron_L4s <- data.table::fread(
  file = file.path(rscdir, "AllenIMNTopMajorRegionRelatedL4s.txt"),
  header = FALSE)$V1
neuronATACSeu <- readRDS(
  file.path(outdir, "neuron_atac_noraw_seurat.rds")
)
imnNeuronAllenSeu <- subset(neuronATACSeu, subset = L4 %in% imn_neuron_L4s)
saveRDS(imnNeuronAllenSeu,
  file = file.path(outdir, "imneuron_atac_noraw_seurat.rds"))
#===================================================================

# * test
outdir <- file.path(here::here(), "19.snap2_integration", "out/test/dp")
dir.create(outdir)

annfnm <- file.path(rscdir, "norawdata_allen", "nn_male_10xv3_ann_noraw.h5ad")
## NOTE: Seurat FindVariableFeatures needs count layer
## and works for Seurat v5
allenS5 <- convertAnn2Seurat5(
  annfnm = annfnm,
  modality = "rna",
  group = "X",
  outdir = file.path(outdir, "nn_10xv3.s5"),
  isLogNorm = TRUE,
  removeCounts = FALSE,
  overwrite =  FALSE,
  addObsMeta = TRUE,
  usePython = "/home/szu/mambaforge/envs/sa2/bin/python"
)

# - this will copy data to the destdir
#   with another directory with the same name is created
## saveRDS(allenS5,
##   file = "s5.rds",
##   destdir = file.path(outdir, "nn_10xv3.s5"))

# - if keep the original one, it will try to save but with
#   error, file are existed.
## saveRDS(allenS5, file = "s5.rds", destdir = outdir)

# in mem
dp.allenS5 <- downsampleSeurat(
  seu = allenS5, groupBy = "cl", minNum = 50, maxNum = 50) |>
  toSeuratInMemory(seu = _, slot = "data", removeCounts = FALSE)
saveRDS(dp.allenS5, file.path(outdir, "dp.allenS5.in.mem.rds"))

dp2.allenS5 <- downsampleSeurat(
  seu = allenS5, groupBy = "cl", minNum = 50, maxNum = 50)
saveRDS(dp2.allenS5, file.path(outdir, "dp.allenS5.not.in.mem.rds"))


# s5 format
dp.allenS5 <- downsampleSeurat(
  seu = allenS5, groupBy = "cl", minNum = 50, maxNum = 50)



# - subset will save all the raw data
#   instead of only subsetting ones
# - under that directory, another directory
#   with the same original dir name will be created.   
## saveRDS(dp.allenS5, file = "dp.allenS5.rds",
  ## destdir = file.path(outdir, "dp.allenS5"))

# * test integration with demo dataset
dp.allenS5$split.id <- sample(
  1:2, size = ncol(dp.allenS5), replace = TRUE)
dplist <- SplitObject(dp.allenS5, split.by = "split.id")
features <- geneList$TF

dp1 <- dplist[[1]]  |>
  ScaleData(object = _, features = features) |>
  RunPCA(object = _, features = features, verbose = TRUE, npcs = 30)
dp1$modality <- "RNA"

dp2 <- dplist[[2]] |>
  ScaleData(object = _, features = features) |>
  RunPCA(object = _, features = features, verbose = TRUE, npcs = 30)
dp2$modality <- "ATAC"

dpMerge <- merge(dp1, dp2,
  add.cell.ids =  NULL,
  merge.data = TRUE, merge.dr = TRUE)

# Integration
# Seurat Integration Reference:
# https://satijalab.org/seurat/archive/v3.0/integration.html

# IntegrateLayers will call method directly
# For each integration
# 1. call FindIntegrationAnchors to get anchors
# 2. then call IntegrateEmbeddings
#    - using anchors to transfer queries embedding to
#      the reference ones, while keep reference unchanged
#    - anchors will be ignroed in the results
#    - reference will be determinied as the largest set if
#      no specified


# ** rpca
# Bug: https://github.com/satijalab/seurat/issues/7720
ints5 <- IntegrateLayers(
  object = dpMerge,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "intgn.rpca",
  verbose = TRUE,
  features = features,
  dims = 1:30,
  k.anchor = 5,
  group.by = "modality",
  scale = FALSE,
  normalization.method = "LogNormalize"
)

# rpca will do not consider scale, will
# run normalziation and PCA again for each
# dataset seperately

## Computing within dataset neighborhoods
## debug: k.neighbor <- max(k.anchor, k.score)
## - k.score = 30, k.anchor = 5
## - k.filter will be set as NA
anchors <- Seurat::FindIntegrationAnchors(
  object.list = list(dp1, dp2),
  anchor.features = features,
  scale = TRUE,
  reduction = "rpca",
  normalization.method = "LogNormalize",
  dims = 1:30,
  k.filter = NA,
  reference = 1,
  verbose = TRUE
)

## Transfer labels use transfer anchors
## rpca works for both in-memory seurat and on-disk s5 BPCell format
tfAnchors <- Seurat::FindTransferAnchors(
  reference = dp1,
  # will not runPCA on ref
  # but runPCA on query needs every time.
  reference.reduction = "pca",
  query = dp2,
  normalization.method = "LogNormalize",
  recompute.residuals = TRUE,
  reduction = "rpca",
  scale = FALSE,
  npcs = 30,
  features = features,
  dims = 1:30,
  k.anchor = 20
)

## Use IntegrateData to get co-embedding

labelstf <- Seurat::TransferData(anchorset = tfAnchors,
  refdata = as.character(dp1$cl),
  query = dp2,
  dims = 1:30,
  weight.reduction = "rpca.ref.l2")

predIds <- labelstf$predicted.id
predScoreMat <- labelstf@assays$prediction.score.id
sum(predIds == as.character(dp2$cl))

## we may can use ths random sampling procedure to choose
## parameters

# ** test CCA
# Running CCA
# Error: Not compatible with requested type: [type=S4; target=double].
# directoly on Seurat v5

tfAnchorCCA <- Seurat::FindTransferAnchors(
  reference = dp1,
  # reference reduction only for 'pcaproject'
  ## reference.reduction = "pca",
  query = dp2,
  normalization.method = "LogNormalize",
  recompute.residuals = TRUE,
  reduction = "cca",
  scale = FALSE,
  npcs = 50,
  dims = 1:30,
  features = features,
  k.anchor = li20
)

labelstf <- Seurat::TransferData(anchorset = tfAnchorCCA,
  refdata = as.character(dp1$cl),
  query = dp2,
  dims = 1:30,
  ## must be cca
  weight.reduction = "cca.l2"
  )

predIds <- labelstf$predicted.id
predScoreMat <- labelstf@assays$prediction.score.id
sum(predIds == as.character(dp2$cl))
