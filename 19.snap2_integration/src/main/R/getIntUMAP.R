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

# * load snakemake configs
intS5fnm <- snakemake@input[[1]]
outfnm <- snakemake@output[[1]]
nPCA <- snakemake@params$nPCA
intMethod <- snakemake@wildcards$m

# * set logger
setupLogging(snakemake@log[[1]])

# * main
logger::log_info("readRDS: ", intS5fnm)
seu <- readRDS(intS5fnm)
logger:log_info("After integration with: ", intMethod)
logger::log_info("run UMAP on reduction: ", paste0("intgn.", intMethod))

seu <- Seurat::RunUMAP(seu,
  reduction = paste0("intgn.", intMethod), dims = 1:nPCA,
  reduction.name = paste0("umap.", intMethod))

logger::log_info("finish UMAP, and save to ", outfnm)
saveRDS(seu, outfnm)
logger::log_info("done")
closeLogging()
