library(Seurat)
options(Seurat.object.assay.version = "v5")

# Not sure if this is a must
library(SeuratObject)
# Not sure if this is a must
library(SeuratData)

# seurat5 version wrappers
# needed for FastMNNIntegration
library(SeuratWrappers)

# [FIXED] Seurat5 use future inside, this is a must
options(future.globals.maxSize = 5e9)

library(stringr)
library(purrr)
library(Matrix)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  Sa2Integration)
import::from(.from = "utils.R", .directory = rdir,
  setupLogging, closeLogging)

# * load snakemake configs
atacS5fnm <- snakemake@input$atacS5
allenS5fnm <- snakemake@input$allenS5
intS5fnm <- snakemake@output[[1]]
logfnm <- snakemake@log[[1]]

nPCA <- snakemake@params$nPCA
fromATAC <- snakemake@params$fromATAC
fromRNA <- snakemake@params$fromRNA

ft <- snakemake@wildcards$ft
kAnchor <- as.integer(snakemake@wildcards$k)
intMethod <- snakemake@wildcards$m

# * set logger
setupLogging(logfnm)

# * load features from Sa2Integration
logger::log_info("load features from Sa2Integration: ", ft)
geneList <- Sa2Integration$getMarkersList()
fts <- geneList[[ft]]

# * main
logger::log_info("readRDS: ", atacS5fnm)
seuATAC <- readRDS(atacS5fnm)
logger::log_info("readRDS: ", allenS5fnm)
seuRNA <- readRDS(allenS5fnm)

## intS5 <- merge(seuATAC, seuRNA,
##   add.cell.ids = NULL, merge.data = TRUE, merge.dr = TRUE)
intS5 <- merge(seuRNA, seuATAC
  add.cell.ids = NULL, merge.data = TRUE, merge.dr = TRUE)

logger::log_info("Run integration with ",
  intMethod, " and k = ", kAnchor)

method <- if (intMethod == "cca") {
  CCAIntegration
} else if (intMethod == "rpca") {
  RPCAIntegration
} else if (intMethod == "mnn") {
  FastMNNIntegration
} else if (intMethod == "harmony") {
  HarmonyIntegration
} else {
  stop("unknown integration method: ", intMethod)
}

intS5 <- IntegrateLayers(
  object = intS5,
  method = method,
  orig.reduction = "pca",
  new.reduction = paste0("intgn.", intMethod),
  verbose = TRUE,
  features = fts,
  scale.layer = "scale.data",
  ## FIXME: if set. will fail.
  ## scale = FALSE,
  normalization.method = "LogNormalize",
  group.by = "modality",
  dims = fromRNA:nPCA,
  k.anchor = kAnchor
)

logger::log_info("finish integration and now save it to: ", intS5fnm)
saveRDS(intS5, intS5fnm)
logger::log_info("done")
closeLogging()
