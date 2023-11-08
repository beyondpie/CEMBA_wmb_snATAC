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

# load my scripts
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  Sa2Integration, cluSumBySa2)
import::from(.from = "ggtheme.R", .directory = rdir,
  dotplotTheme)

import::from(.from = "integration.R", .directory = rdir,
  TransferLabelSum, getPredictLabel,
  getPredictScoreMat, getMetaCol.Seurat,
  prepareDotPlot4TransferLabel)

import::from(.from = "utils.R", .directory = rdir,
  toNamedArray.1dtable, to3.matrix,
  getMaxColScoreWithName)

mytheme <- dotplotTheme(legend.pos = "right")

# * main
# load tfSum
fromdir <- file.path(here::here(), "19.snap2_integration",
  "out/tf_tscc/neuron")

kAnchor <- 30
feaName <- "k8"
allenTech <- "allen"
method <- "cca"

tfSum <- TransferLabelSum$new(
  kAnchor = kAnchor,
  feaName = feaName,
  dsNum = 50,
  allenTech = allenTech,
  group = "neuron",
  method = method
)
tfSum$loadtfLabel(fromdir = fromdir)
refCol <- "subclass_id"
queryCol <- "L4"
tfVoteMat <- tfSum$getAvgVoteMat.queryGroup(tfSum,
  refCol = refCol,
  queryCol = queryCol)
mat.dotplot <- prepareDotPlot4TransferLabel(
  tfmat = tfVoteMat, ignoreEmptyRef = TRUE,
  names = c("Allen", "ATAC", "score"))
maxScore <- getMaxColScoreWithName(tfVoteMat)
query2ref_subclassid <- data.frame(
  query = colnames(tfVoteMat),
  ref = names(maxScore),
  row.names = colnames(tfVoteMat)
)

# load Allen meta
allenAnnotMeta <- Sa2Integration$loadAllenAnnot()

refCol <- "cl"
queryCol <- "L4"
tfVoteMat <- tfSum$getAvgVoteMat.queryGroup(tfSum,
  refCol = refCol,
  queryCol = queryCol)
mat.dotplot <- prepareDotPlot4TransferLabel(
  tfmat = tfVoteMat, ignoreEmptyRef = TRUE,
  names = c("Allen", "ATAC", "score"))

maxScore <- getMaxColScoreWithName(tfVoteMat)
query2ref <- data.frame(
  query = colnames(tfVoteMat),
  ref = names(maxScore),
  row.names = colnames(tfVoteMat)
)

length(unique(allenAnnotMeta[unique(query2ref$ref), "subclass_id"]))
length(unique(allenAnnotMeta[unique(query2ref$ref), "supertype_id"]))

query2ref$subclass_id <- allenAnnotMeta[query2ref$ref, "subclass_id"]
query2ref$subclass_id_tf <- as.integer(
  query2ref_subclassid[query2ref$query, "ref"])
a <- with(query2ref, subclass_id == subclass_id_tf)

length(unique(query2ref[!a, "subclass_id_tf"]))

sum(allenAnnotMeta[unique(query2ref$ref), "v3.size"])
sum(allenAnnotMeta[
  allenAnnotMeta$subclass_id %in% query2ref$subclass_id, "v3.size"])
sum(allenAnnotMeta[rownames(tfVoteMat), "v3.size"]) 


# load snapatac2 clustering info
barcode2L4 <- cluSumBySa2$loadbarcode2L4()
L4size <- table(barcode2L4$L4)
barcode2L3 <- cluSumBySa2$loadbarcode2L3()
L3size <- table(barcode2L3$L3)
