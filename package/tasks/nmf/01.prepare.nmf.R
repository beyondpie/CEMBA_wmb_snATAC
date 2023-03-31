# Prepare data for nmf.
# TODO: This script can be generalised.
# Currently, I use it case by case.
library(data.table)
library(hdf5r)

packdir <- file.path(here::here(), "package", "R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env)

# * configs
outdir <- "data"
outh5 <- file.path(outdir, "cpm.cbyp.Intv2.h5")
outPeakCoordFile <- file.path(outdir, "peaks.Intv2.txt")
outClusterFile <- file.path(outdir, "clusters.Intv2.txt")

# * load atac cpm.
cpm.scbyp <- readRDS(cembav2env$subclassPmatCPMIntv2File)

## # * load cCREs from posivite pdc
## cCREs.ppdc <- data.table::fread(
##   "../out/AllenAnnotConcat/mba.whole.AllenAnnotConcat.pearson.pos.pdc.CREs",
##   header = FALSE, data.table = FALSE
## )$V1

# * main
## cpm.ppdc <- cpm.pdc[cCREs.ppdc, ]
## cpm.pbysc <- cpm.ppdc

## set cap for too high values
upValue <- quantile(cpm.scbyp, 0.9999)
cpm.capped <- cpm.scbyp
cpm.capped[cpm.scbyp > upValue] <- upValue

## change to cluster by peak mat
## for later saving to hdf5
peaks <- colnames(cpm.capped)
clusters <- rownames(cpm.capped)

# * save mat to hdf5 format for python handling.
conn <- hdf5r::H5File$new(outh5, mode = "w")
data.grp <- conn$create_group("X")
# NOTE: hdf5r will transpose the mat
# https://github.com/hhoeflin/hdf5r/issues/81
data.grp[["mat"]] <- cpm.capped
# colnames corresponds to cpm.capped
data.grp[["colnames"]] <- peaks
# rownames corresponds to cpm.capped
data.grp[["rownames"]] <- clusters
conn$close_all()

write.table(peaks, file = outPeakCoordFile, quote = FALSE,
  row.names = FALSE, col.names = FALSE)
write.table(clusters, file = outClusterFile, quote = FALSE,
  row.names = FALSE, col.names = FALSE)
