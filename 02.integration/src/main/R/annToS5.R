library(BPCells)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(stringr)
library(purrr)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "utils.R", .directory = rdir,
  setupLogging, closeLogging)


## library(future)
## plan(multicore, workers = 2)
## options(future.globals.maxSize = 2e9)
## options(future.rng.onMisuse = "ignore")


# * load snakemake configs
annfnm <- snakemake@input[[1]]
outfnm <- snakemake@output[[1]]
logfnm <- snakemake@log[[1]]
modality <- snakemake@params[["modality"]]
## py <- snakemake@params[["py"]]

# * set logger
setupLogging(logfnm)

# * set python
## library(reticulate)
## reticulate::use_python(py)
## ad <- reticulate::import("andnata", convert = FALSE)

# * function
convertAnn2Seurat5 <- function(annfnm,
                               modality,
                               group = "X",
                               outdir,
                               overwrite = TRUE,
                               assay = "RNA",
                               isLogNorm = TRUE) {
  xlognorm <- BPCells::open_matrix_anndata_hdf5(
    path = annfnm, group = group)
  BPCells::write_matrix_dir(mat = xlognorm, dir = outdir, overwrite = overwrite)
  d <- BPCells::open_matrix_dir(outdir)
  s5 <- Seurat::CreateSeuratObject(counts = d, assay = assay)
  if (isLogNorm) {
    s5 <- Seurat::SetAssayData(
      object = s5, slot = "data", new.data = d)
  }
  s5$modality <- modality
  return(s5)
}


# * main
outdir <- dirname(outfnm)
outfprefix <- sub("\\.[^.]+$", "", basename(outfnm))
logger::log_info("output to: ", outdir)
logger::log_info("out file prefix: ", outfprefix)

logger::log_info("to Seuratv5: ", annfnm)
s5 <- convertAnn2Seurat5(
  annfnm = annfnm, modality = modality,
  outdir = file.path(outdir, outfprefix), isLogNorm = TRUE)
logger::log_info(
  "NOTE: [modality] column added to Seuratv5 with value: ", modality)
logger::log_info("NOTE: obs meta data is ignored in Seuratv5.")

logger::log_info("save seuratv5 to: ", outfnm)
saveRDS(s5, outfnm)
logger::log_info("to Seuratv5 done.")
closeLogging()
