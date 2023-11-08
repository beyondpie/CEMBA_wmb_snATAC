# save all the critical meta data into one file
# can be shard with the community
# and used for figures generation
library(tidyverse)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratObject)

library(stringr)
library(purrr)
# support R new feature in 4.2
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

library(R6)
library(future)
options(future.globals.maxSize = 10e9)
library(future.apply)
plan(multicore, workers = 4)

library(ggplot2)

# * sample-level data
# [DATA] organize into one data.frame
# add replicate info, add sample-level stats info
# rep: early, later

# REF:
# - meta/dissect2time.csv
# - meta/BrainRegion.Metadata.txt
# - fisrt submission meta:
#   - supple.02.annotation.all.in.one.meta/mba.whole.cell.meta.v8.1.rds

# * single-cell level
# [DATA] Single-cell level can be saved into one big dataframe
# clustering results:
# REF:
# /projects/ps-renlab2/szu/projects/CEMBA2/17.snapatac2/post_script/17.snapatac2/post_script/sa2.barcode2L4.csv

# ** from snapatac2 sample h5ad, read tsse, #fragment, dlt score
#    and all the barcodes passing QC
#    this is done in python
# REF:
# Before Doublet: sample-level, h5ad, snapatac2
# /oasis/tscc/scratch/szu/projects/CEMBA2/17.snapatac2/sa2_qc_dlt
# After Doublet: sample-level, h5ad, snapatac2
# /projects/ps-renlab2/szu/projects/CEMBA2/17.snapatac2/sa2_qc_dlt/rm_dlt
# Code Ref:
## reticulate::use_python(usePython)
## ad <- reticulate::import("anndata", convert = FALSE)
## allen_ann <- ad$read_h5ad(filename = annfnm, backed = 'r')
## obs_meta <- allen_ann$obs
## obsMeta <- reticulate::py_to_r(obs_meta)
## attr(obsMeta, "pandas.index") <- NULL
## if ("cl" %in% colnames(obsMeta)) {
##   cls <- obsMeta$cl
##   if (is.factor(cls)) {
##     message("Hard code 'cl' column from factor to int.")
##     clInt <- as.integer(levels(cls)[cls])
##     obsMeta$cl <- clInt
##   }
## } # end of transform cl to Int


# ** based on original meta, add rough annot for our barcodes
# REF:
# 19.snap2_integration/src/main/resource/sa2_dlt2_L3toSa1Annot.rough.csv

# * organize Allen's latest meta about annotation
# [DATA] Allen's can be saved into one data.frame
# REF:
# - CEMBA2/19.snap2_integration/src/main/resource/AIT21_annotation_freeze_081523.tsv
# - AIT21_k8_markers.txt

# * cluster level
# [DATA] all the cluster level info can be saved into oen data.frame

# ** L4 to Allen's cl, supertype, subclass, class
# add a tag to this version: v1
# REF:
# /projects/ps-renlab2/szu/projects/CEMBA2/19.snap2_integration/out/tf_lightsum
# - neuron_atac.allen_cca_k8_50.lightsum.rds
# - nn_atac.allen_cca_k8_49.lightsum.rds

# ** map major and sub region information to L4-class

# ** get allen subclass region info and mapped to CEMBA
# [DATA] extra file for checking the integration
# mark with different versions
# REF:
# meta/allen.region.to.main.region.txt

# * check result from Kangli
atacMeta <- readRDS(
  file.path(here::here(),
    "supple.02.annotation.all.in.one.meta",
    "mba.whole.cell.meta.v9.2.rds"))

# check balance
l4 <- "6-5-2-3"
l4 <- "10-3-1-2"
table(with(atacMeta, biorep[L4 == l4]))

# check subclass annotation
# compare with my previous one.
a <- with(atacMeta, table(AllenAnnot.rough, subclass_label.allen))
b <- a / rowSums(a)

b <- as.data.frame(b)
c <- b[!(b$Freq < 0.0001), ]

(
  p <- ggplot(data = c, aes(x = subclass_label.allen, y = AllenAnnot.rough)) +
  geom_point(aes(color = Freq)) +
  scale_color_gradient(low = "white", high = "red", na.value = "white") +
  theme(
    axis.text.x = element_text(color = "black", size = 4.5, angle = 90),
    axis.text.y = element_text(color = "black", size = 4.5)
  )
)

# * get subclass 2 size csv.
atacMeta <- readRDS(
  file.path(here::here(),
    "supple.02.annotation.all.in.one.meta",
    "mba.whole.cell.meta.v9.2.rds"))

tmp <- table(atacMeta$subclass_id_label.allen) |>
  as.data.frame() |>
  mutate(.data = _, subclass = gsub(" +", "_", as.character(Var1))) 

write.table(tmp[, c("subclass", "Freq")],
  file = file.path(here::here(), "19.snap2_integration",
    "src/main/resource/atac.subclass2size.v1.csv"),
  sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)

# * update subclass annotation after manually checking
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env, Sa2Integration)
atacMeta <- readRDS(cembav2env$sa2metaFile)

rscdir <- file.path(projdir, "19.snap2_integration",
  "src/main/resource")
# subclass_id_v2 based on 230928 version
# now upgrade to v3
scManualCheckfnm <- file.path(
  rscdir, "atac_top3_subclass_checkRegion_230930.tsv")
scchecked <- data.table::fread(
  file = scManualCheckfnm,
  sep = "\t",
  header = TRUE,
  data.table = FALSE)
rownames(scchecked) <- scchecked$L4
atacMeta$subclass_id_v3 <- scchecked[atacMeta$L4, "subclass_id"]
allenAnnot <- Sa2Integration$loadAllenAnnot()
sc_id2nm <- unique(allenAnnot[,
  c("subclass_id", "subclass_label", "subclass_id_label")])

rownames(sc_id2nm) <- as.character(sc_id2nm$subclass_id)
atacMeta$subclass_id_label_v3 <- sc_id2nm[
  as.character(atacMeta$subclass_id_v3), "subclass_id_label"]

atacMeta$subclass_label_v3 <- sc_id2nm[
  as.character(atacMeta$subclass_id_v3), "subclass_label"]
saveRDS(atacMeta, file.path(
  projdir, "supple.02.annotation.all.in.one.meta",
  "mba.whole.cell.meta.v9.6.rds"))

# * check OGC
