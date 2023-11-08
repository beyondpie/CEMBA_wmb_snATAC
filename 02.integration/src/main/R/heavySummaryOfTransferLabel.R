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

library(Matrix)

projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  Sa2Integration, cluSumBySa2)
import::from(.from = "integration.R", .directory = rdir,
  TransferLabelSum, getPredictLabel,
  getPredictScoreMat, getMetaCol.Seurat,
  prepareDotPlot4TransferLabel)
import::from(.from = "utils.R", .directory = rdir,
  toNamedArray.1dtable, to3.matrix,
  getMaxColScoreWithName,
  getTopkColScoreWithName)

# * main
allenMeta <- Sa2Integration$loadAllenAnnot()
cl2subclass <- allenMeta[, c("cl", "subclass_id")]
rownames(cl2subclass) <- cl2subclass$cl
sc2region <- Sa2Integration$loadAllenSubclass2Region()

workdir <- file.path(projdir, "19.snap2_integration",
  "out/tf_tscc")
outdir <- file.path(projdir, "19.snap2_integration",
  "out/topktf")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# get top 3 matched ones.
topk <- 3L
refCol <- "cl"
getSubclasstfScoreByclId <- function(scMat, top3cl, cl2subclass) {
 lapply(seq_along(top3cl), \(i) {
    l4 <- names(top3cl)[i]
    top3 <- top3cl[[i]]
    subclasses <- cl2subclass[names(top3), "subclass_id"] |>
      unique()
    scMat[as.character(subclasses), l4] |>
      setNames(object = _, nm = subclasses)
  }) |>
    setNames(object = _, nm = names(top3cl))
}
feaName <- "k8"
allenTech <- "allen"
method <- "cca"

## neuron
kAnchor <- 50L
group <- "neuron"
dsNum <- 50
tffnm <- file.path(workdir, group,
  str_glue("{group}_atac.{allenTech}_{method}_{feaName}_{kAnchor}.tf.anchor.rds"))
# outfnm <- file.path(outdir, str_glue("{group}_top{topk}.csv"))

neurontfSum <- TransferLabelSum$new(
  kAnchor = kAnchor,
  feaName = feaName,
  dsNum = dsNum,
  allenTech = allenTech,
  group = group,
  method = method
)
# cost lots of memory (70G)
neurontfSum$loadtfLabel(fromdir = file.path(workdir, group))

neuron.scMat <- neurontfSum$getAvgVoteMat4SubclassByRefclId(
  tfSum = neurontfSum, queryCol = "L4", cl2subclass = cl2subclass)
saveRDS(neuron.scMat,
  file.path(outdir, str_glue("{group}_subclass2L4_tfscore.rds")))

neuron.clMat <- neurontfSum$getAvgVoteMat.queryGroup(
  tfSum = neurontfSum, refCol = "cl", queryCol = "L4")
saveRDS(neuron.clMat,
  file.path(outdir, str_glue("{group}_cl2L4_tfscore.rds")))

neuron.top3cl <- getTopkColScoreWithName(
  neuron.clMat, topk = 3, decreasing = TRUE)
neuron.top3sc <- getSubclasstfScoreByclId(
  neuron.scMat, neuron.top3cl, cl2subclass)
saveRDS(neuron.top3sc,
  file.path(outdir, str_glue("{group}_top{topk}.rds")))
rm(neurontfSum)

## rnames <- rownames(neuron.tfVoteMat)
## tmp <- neuron.tfVoteMat[,1]
## tmp[order(tmp, decreasing = TRUE)][1:3]
## rnames[order(tmp, decreasing = TRUE)][1:3]

# * nn
kAnchor <- 49L
group <- "nn"
dsNum <- 500
tffnm <- file.path(workdir, group,
  str_glue("{group}_atac.{allenTech}_{method}_{feaName}_{kAnchor}.tf.anchor.rds"))
# outfnm <- file.path(outdir, str_glue("{group}_top{topk}.csv"))

nntfSum <- TransferLabelSum$new(
  kAnchor = kAnchor,
  feaName = feaName,
  dsNum = dsNum,
  allenTech = allenTech,
  group = group,
  method = method
)
nntfSum$loadtfLabel(fromdir = file.path(workdir, group))
nn.scMat <- nntfSum$getAvgVoteMat4SubclassByRefclId(
  tfSum = nntfSum, queryCol = "L4", cl2subclass = cl2subclass)
saveRDS(nn.scMat,
  file.path(outdir, str_glue("{group}_subclass2L4_tfscore.rds")))

nn.clMat <- nntfSum$getAvgVoteMat.queryGroup(
  tfSum = nntfSum, refCol = "cl", queryCol = "L4")
saveRDS(nn.clMat, file.path(outdir, str_glue("{group}_cl2L4_tfscore.rds")))
neuron.top3cl <- getTopkColScoreWithName(
  nn.clMat, topk = 3, decreasing = TRUE)
nn.top3sc <- getSubclasstfScoreByclId(
  nn.scMat, neuron.top3cl, cl2subclass)
saveRDS(nn.top3sc,
  file.path(outdir, str_glue("{group}_top{topk}.rds")))
rm(nntfSum)


# * imn
## imneuron
kAnchor <- 50L
group <- "imneuron"
dsNum <- 50
tffnm <- file.path(workdir, group,
  str_glue("{group}_atac.dp50allen_{method}_{feaName}_{kAnchor}.tf.anchor.rds"))
imntfSum <- TransferLabelSum$new(
  kAnchor = kAnchor,
  feaName = feaName,
  dsNum = dsNum,
  allenTech = "dp50allen",
  group = group,
  method = method
)
imntfSum$loadtfLabel(fromdir = file.path(workdir, group))
imn.clMat <- imntfSum$getAvgVoteMat.queryGroup(
  tfSum = imntfSum, refCol = "cl", queryCol = "L4")
saveRDS(imn.clMat,
  file.path(outdir, str_glue("{group}_cl2L4_tfscore.rds")))
imn.scMat <- imntfSum$getAvgVoteMat4SubclassByRefclId(
  tfSum = imntfSum, queryCol = "L4", cl2subclass = cl2subclass)
saveRDS(imn.scMat,
  file.path(outdir, str_glue("{group}_subclass2L4_tfscore.rds")))

# * merge cl to L4 mat
nn.clMat <- readRDS(file.path(
  outdir,
  "nn_cl2L4_tfscore.rds"
))
neuron.clMat <- readRDS(file.path(
  outdir,
  "neuron_cl2L4_tfscore.rds"
))
imn.clMat <- readRDS(file.path(
  outdir,
  "imneuron_cl2L4_tfscore.rds"
))

neuron2nn.clMat <- matrix(0.0,
  nrow = nrow(neuron.clMat), ncol = ncol(nn.clMat))
colnames(neuron2nn.clMat) <- colnames(nn.clMat)
rownames(neuron2nn.clMat) <- rownames(neuron.clMat)

neuron.all.clMat <- cbind(neuron.clMat, neuron2nn.clMat)
rownames(neuron.all.clMat) <- rownames(neuron.clMat)
colnames(neuron.all.clMat) <- c(
  colnames(neuron.clMat),
  colnames(neuron2nn.clMat)
)

nn2neuron.clMat <- matrix(0.0,
  nrow = nrow(nn.clMat), ncol = ncol(neuron.clMat))
colnames(nn2neuron.clMat) <- colnames(neuron.clMat)
rownames(nn2neuron.clMat) <- rownames(nn.clMat)

nn.all.clMat <- cbind(nn2neuron.clMat, nn.clMat)
rownames(nn.all.clMat) <- rownames(nn.clMat)
colnames(nn.all.clMat) <- c(
  colnames(neuron.clMat),
  colnames(nn.clMat)
)

all.clMat <- rbind(neuron.all.clMat, nn.all.clMat)
rownames(all.clMat) <- c(
  rownames(neuron.all.clMat),
  rownames(nn.all.clMat)
)
colnames(all.clMat) <- colnames(nn.all.clMat)

# update imn to neuron
cl_imn <- allenMeta$cl[grep("IMN", allenMeta$subclass_label)]
cl_imn_str <- as.character(cl_imn)
all(cl_imn_str %in% rownames(imn.clMat))
imn2L4.clMat <- imn.clMat[cl_imn_str, ]

all.clMat[rownames(imn2L4.clMat), colnames(imn2L4.clMat)] <- imn2L4.clMat
saveRDS(all.clMat, file.path(outdir, "sa2.all.cl2L4.mat.rds"))
write.table(all.clMat, file.path(outdir, "sa2.all.cl2L4.transferLabelScore.csv"),
  sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)


# * merge sc to L4 mat
## NOTE: we use the variable name before but load sc
nn.clMat <- readRDS(file.path(
  outdir,
  "nn_subclass2L4_tfscore.rds"
))
neuron.clMat <- readRDS(file.path(
  outdir,
  "neuron_subclass2L4_tfscore.rds"
))
imn.clMat <- readRDS(file.path(
  outdir,
  "imneuron_subclass2L4_tfscore.rds"
))

neuron2nn.clMat <- matrix(0.0,
  nrow = nrow(neuron.clMat), ncol = ncol(nn.clMat))
colnames(neuron2nn.clMat) <- colnames(nn.clMat)
rownames(neuron2nn.clMat) <- rownames(neuron.clMat)

neuron.all.clMat <- cbind(neuron.clMat, neuron2nn.clMat)
rownames(neuron.all.clMat) <- rownames(neuron.clMat)
colnames(neuron.all.clMat) <- c(
  colnames(neuron.clMat),
  colnames(neuron2nn.clMat)
)

nn2neuron.clMat <- matrix(0.0,
  nrow = nrow(nn.clMat), ncol = ncol(neuron.clMat))
colnames(nn2neuron.clMat) <- colnames(neuron.clMat)
rownames(nn2neuron.clMat) <- rownames(nn.clMat)

nn.all.clMat <- cbind(nn2neuron.clMat, nn.clMat)
rownames(nn.all.clMat) <- rownames(nn.clMat)
colnames(nn.all.clMat) <- c(
  colnames(neuron.clMat),
  colnames(nn.clMat)
)

all.clMat <- rbind(neuron.all.clMat, nn.all.clMat)
rownames(all.clMat) <- c(
  rownames(neuron.all.clMat),
  rownames(nn.all.clMat)
)
colnames(all.clMat) <- colnames(nn.all.clMat)

# update imn to neuron
cl_imn <- allenMeta$subclass_id[grep("IMN", allenMeta$subclass_label)]
cl_imn_str <- as.character(cl_imn)
all(cl_imn_str %in% rownames(imn.clMat))
imn2L4.clMat <- imn.clMat[cl_imn_str, ]

all.clMat[rownames(imn2L4.clMat), colnames(imn2L4.clMat)] <- imn2L4.clMat
saveRDS(all.clMat, file.path(outdir, "sa2.all.subclass2L4.mat.rds"))
write.table(all.clMat, file.path(outdir,
  "sa2.all.subclass2L4.transferLabelScore.csv"),
  sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)

# * save cl to subclass id to subclass label
cl2sc2label <- unique(allenMeta[ , c("cl", "subclass_id", "subclass_id_label")])
write.table(cl2sc2label,
  file.path(outdir, "Allen.scRNAseq.cl_subclassid_subclasslabel.csv"),
  sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
# remove LQ manually
