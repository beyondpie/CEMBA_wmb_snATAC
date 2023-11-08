library(hdf5r)
library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(future)
library(data.table)

projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env, Sa2PeakCalling)
import::from(.from = "utils.R", .directory = rdir,
  fastread.csv)

workdir <- file.path(projdir, "20.nmf")
peaksumdir <- file.path(projdir, "supple.07.peakcalling.allinone")

cpm_pbysc <- fastread.csv(
  fnm = file.path(projdir, "18.snap2_peakcalling",
    "out/scfilter",
    "cpm_peakBysubclass.csv"))
cpm_pbysc_mat <- as.matrix(cpm_pbysc)
upValue <- quantile(cpm_pbysc_mat, 0.9999)
cpm_capped <- cpm_pbysc_mat
cpm_capped[cpm_capped > upValue] <- upValue
# to subclass by peak
cpm_capped <- t(cpm_capped)
peaks <- colnames(cpm_capped)
clusters <- rownames(cpm_capped)

# * save the pmat with all the peak into a hdf5 format
conn <- hdf5r::H5File$new(
  file.path(workdir, "data/sa2.allpeak.cbyp.mat.h5"), mode = "w")
data.grp <- conn$create_group("X")
# NOTE: hdf5r will transpose the mat
# https://github.com/hhoeflin/hdf5r/issues/81
data.grp[["mat"]] <- cpm_capped
# colnames corresponds to cpm.capped
data.grp[["colnames"]] <- peaks
# rownames corresponds to cpm.capped
data.grp[["rownames"]] <- clusters
conn$close_all()

write.table(peaks,
  file = file.path(workdir, "data/sa2.final.peak.nms.txt"),
  quote = FALSE,
  row.names = FALSE, col.names = FALSE)
write.table(clusters,
  file = file.path(workdir, "data/sa2.allpeak.subclass.nms.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

# * set order of the clusters
sa2SubclassNames <- data.table::fread(
  file = file.path(projdir, "meta", "sa2.subclass.names.map.csv"),
  header = TRUE, sep = ",", data.table = FALSE)
rownames(sa2SubclassNames) <- sa2SubclassNames$subclass_label_peak
subclass_ids <- sa2SubclassNames[clusters, "subclass_id"]
sc_ordered <- clusters[order(subclass_ids, decreasing = FALSE)]
write.table(sc_ordered,
  file = file.path(workdir, "data/sa2.subclass.srt.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

# * prepare the hdf5 file for ppdc
peak.ppdc <- data.table::fread(
  file = cembav2env$sa2.pos.pdc.peak.nm,
  header = FALSE, data.table = FALSE
)$V1
cpm_ppdc <- cpm_capped[ , peak.ppdc]
conn_ppdc <- hdf5r::H5File$new(
  file.path(workdir, "data/sa2.ppdc.cbyp.mat.h5"), mode = "w")
data_ppdc <- conn_ppdc$create_group("X")
data_ppdc[["mat"]] <- cpm_ppdc
data_ppdc[["colnames"]] <- colnames(cpm_ppdc)
data_ppdc[["rownames"]] <- rownames(cpm_ppdc)
conn_ppdc$close_all()
write.table(colnames(cpm_ppdc),
  file = file.path(workdir, "data/sa2.ppdc.peak.nms.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(cpm_ppdc),
  file = file.path(workdir, "data/sa2.ppdc.subclass.nms.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)


# * prepare the hdf5 file for nonovlp peaks
peak.novlpDHS.df <- Sa2PeakCalling$readBed4File(
  cembav2env$sa2.novlp.rDHS.peak.nm)
cpm_novlpDHS <- cpm_capped[ ,peak.novlpDHS.df$name]
conn_novlpDHS <- hdf5r::H5File$new(
  file.path(workdir, "data/sa2.novlpDHS.cbyp.mat.h5"),
  mode = "w")
data_novlpDHS <- conn_novlpDHS$create_group("X")
data_novlpDHS[["mat"]] <- cpm_novlpDHS
data_novlpDHS[["colnames"]] <- colnames(cpm_novlpDHS)
data_novlpDHS[["rownames"]] <- rownames(cpm_novlpDHS)
conn_novlpDHS$close_all()
write.table(colnames(cpm_novlpDHS),
  file = file.path(workdir, "data/sa2.novlpDHS.peak.nms.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(cpm_novlpDHS),
  file = file.path(workdir, "data/sa2.novlpDHS.subclass.nms.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)
