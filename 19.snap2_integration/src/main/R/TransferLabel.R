library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratObject)
options(future.globals.maxSize = 5e9)
library(stringr)
library(purrr)
# support R new feature in 4.2
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
# load my scripts
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "utils.R", .directory = rdir,
  setupLogging, closeLogging)
import::from(.from = "cembav2env.R", .directory = rdir,
  Sa2Integration)
import::from(.from = "integration.R", .directory = rdir,
  convertAnn2Seurat5, get.downsample.fun, downsampleSeurat,
  toSeuratInMemory,
  isOnDiskMat.Seurat, calVarOfFea.Seurat,
  getvf.Seurat,
  setVariableFeatures)

# * snakemake configs
atacseuf <- snakemake@input[["atacSeu"]]
allenseuf <- snakemake@input[["allenSeu"]]
logfnm <- snakemake@log[[1]]

tflabelfnm <- snakemake@output[["transferLabel"]]
tfanchorfnm <- snakemake@output[["anchor"]]
outdir <- dirname(tfanchorfnm)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# data preprocess params
dsSeuart <- as.integer(snakemake@params[["dsSeurat"]])

atac_dsMin <- as.integer(snakemake@params[["atac_dsMin"]])
atac_dsMax <- as.integer(snakemake@params[["atac_dsMax"]])

allen_dsMin <- as.integer(snakemake@params[["allen_dsMin"]])
allen_dsMax <- as.integer(snakemake@params[["allen_dsMax"]])

atacGroupBy <- snakemake@params[["atacGroupBy"]]
allenGroupBy <- snakemake@params[["allenGroupBy"]]

# for get features
useAllenAsRef <- as.integer(snakemake@params[["useAllenAsRef"]])
nfea <- as.integer(snakemake@params[["nFeatures"]])
vfon <- snakemake@params[["vfon"]]

# transfer data params
nPCA <- as.integer(snakemake@params[["nPCA"]])
tfmethod <- snakemake@wildcards[["m"]]
kAnchor <- as.integer(snakemake@wildcards[["k"]])
feaGroup <- snakemake@wildcards[["ft"]]

# * set logger
setupLogging(logfnm)

# * main
# * load seurat data
message("load atac Seurat: ", atacseuf)
atacSeu <- readRDS(atacseuf)

message("load allen Seurat: ", allenseuf)
allenSeu <- readRDS(allenseuf)

if (dsSeuart > 0) {
  message("require: downsample Seurat.")
  message(str_glue("atac: {atac_dsMin} - {atac_dsMax}."))
  atacSeu <- downsampleSeurat(seu = atacSeu, groupBy = atacGroupBy,
    minNum = atac_dsMin, maxNum = atac_dsMax)
  
  message(str_glue("allen: {allen_dsMin} - {allen_dsMax}."))
  allenSeu <- downsampleSeurat(seu = allenSeu, groupBy = allenGroupBy,
    minNum = allen_dsMin, maxNum = allen_dsMax)
}

message("log: cells from atac seq: ", dim(atacSeu)[2])
message("log: cells from allen seq: ", dim(allenSeu)[2])

if (tfmethod == "cca") {
  message("require: cca needs in memory data.")
  atacSeu <- toSeuratInMemory(atacSeu,
    slot = "data", removeCounts = FALSE)
  allenSeu <- toSeuratInMemory(allenSeu,
    slot = "data", removeCounts = FALSE)
}

# * set reference and query dataset
ref <- allenSeu
query <- atacSeu
if (useAllenAsRef < 1) {
  message("NOTE: use ATAC as reference for variable features.")
  ref <- atacSeu
  query <- allenSeu
}


# * set features
message("prepare features from: ", feaGroup)
geneList <- Sa2Integration$getMarkersList()
prefeatures <- NULL
if (feaGroup %in% names(geneList)) {
  prefeatures <- geneList[[feaGroup]]
}
features <- setVariableFeatures(
  ref = ref,
  query = query,
  features = prefeatures,
  onlayer = "counts",
  nfeatures = nfea,
  on = vfon,
  eps = 0.0001
)

# * run find anchors and save file
message("Run: FindTransferAnchors with method: ", tfmethod)
if (file.exists(tfanchorfnm)) {
  message("Warning: ", tfanchorfnm, " exist, and skip FindTransferAnchors.")
  anchors <- readRDS(tfanchorfnm)
} else {
  anchors <- Seurat::FindTransferAnchors(
    reference = ref,
    query = query,
    scale = TRUE,
    reduction = tfmethod,
    normalization.method = "LogNormalize",
    npcs = nPCA,
    dims = seq_len(nPCA),
    l2.norm = TRUE,
    k.anchor = kAnchor,
    features = features,
    verbose = TRUE,
    )
  message("Save: save anchors to ", tfanchorfnm)
  saveRDS(anchors, tfanchorfnm)
}

# * run transferLabel and save file

## wr <- if (tfmethod == "cca") {
##   "cca.l2"
## } else {
##   "rpca.ref.l2"
## }

wr <- if (tfmethod == "cca") {
  "cca"
} else {
  "pcaproject"
}

message("weight reduction for label transfer: ", wr)

refLabelCols <- c("subclass_id", "supertype_id", "cl")

tfLabels <- purrr::map(
  refLabelCols,
  \(x) {
    message("transfer labels from ", x)
    refdata <- as.character(ref[[x]][,1])
    names(refdata) <- colnames(ref)
    Seurat::TransferData(
      anchorset = anchors,
      refdata = refdata,
      query = query,
      dims = seq_len(nPCA),
      weight.reduction = wr,
      l2.norm = FALSE
    )
  }
)
names(tfLabels) <- refLabelCols

message("Save tfLabels to: ", tflabelfnm)
saveRDS(tfLabels, tflabelfnm)

# * clean logging
closeLogging()
message("Run: finished. Good luck!")
