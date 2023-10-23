library(Seurat)
options(Seurat.object.assay.version = "v5")
library(stringr)
library(purrr)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  Sa2Integration)
import::from(.from = "utils.R", .directory = rdir,
  setupLogging, closeLogging)
import::from(.from = "integration.R", .directory = rdir,
  get.downsample.fun, downsampleSeurat)


# * load snakemake configs
infnm <- snakemake@input[[1]]
outfnm <- snakemake@output[[1]]
ft <- snakemake@wildcards$ft
logfnm <- snakemake@log[[1]]
nPCA <- snakemake@params$nPCA

# * set logger
setupLogging(logfnm)

# * load features from Sa2Integration
logger::log_info("load features from Sa2Integration: ", ft)
geneList <- Sa2Integration$getMarkersList()
fts <- geneList[[ft]]

# * main
logger::log_info("readRDS: ", infnm)
seu <- readRDS(infnm)
logger::log_info("ScaleData on feature: ", ft)
Seurat::VariableFeatures(seu) <- fts
seu <- Seurat::ScaleData(seu, features = fts)

logger::log_info("runPCA")
seu <- Seurat::RunPCA(
  seu, features = fts, verbose = TRUE,
  npcs = nPCA)

logger::log_info("writeRDS: ", outfnm)
saveRDS(seu, outfnm)
logger::log_info("done")
closeLogging()
