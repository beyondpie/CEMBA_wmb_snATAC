library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratObject)

library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(future)
options(future.globals.maxSize = 10e9)
library(future.apply)
plan(multicore, workers = 4)
library(Matrix)

projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  Sa2Integration, cluSumBySa2, cembav2env)
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
allenTech <- "dp50allen"
method <- "cca"

## imneuron
kAnchor <- 50L
group <- "imneuron"
dsNum <- 50
tffnm <- file.path(workdir, group,
  str_glue("{group}_atac.{allenTech}_{method}_{feaName}_{kAnchor}.tf.anchor.rds"))
# outfnm <- file.path(outdir, str_glue("{group}_top{topk}.csv"))

imntfSum <- TransferLabelSum$new(
  kAnchor = kAnchor,
  feaName = feaName,
  dsNum = dsNum,
  allenTech = allenTech,
  group = group,
  method = method
)
# cost lots of memory
imntfSum$loadtfLabel(fromdir = file.path(workdir, group))
imn.scMat <- imntfSum$getAvgVoteMat4SubclassByRefclId(
  tfSum = imntfSum, queryCol = "L4", cl2subclass = cl2subclass)
saveRDS(imn.scMat,
  file.path(outdir, str_glue("{group}_subclass2L4_tfscore.rds")))

imn.clMat <- imntfSum$getAvgVoteMat.queryGroup(
  tfSum = imntfSum, refCol = "cl", queryCol = "L4")
imn.top3cl <- getTopkColScoreWithName(
  imn.clMat, topk = 3, decreasing = TRUE)
imn.top3sc <- getSubclasstfScoreByclId(
  imn.scMat, imn.top3cl, cl2subclass)
saveRDS(imn.top3sc,
  file.path(outdir, str_glue("{group}_top{topk}_subclass.rds")))
saveRDS(imn.top3cl,
  file.path(outdir, str_glue("{group}_top{topk}_cl.rds")))
rm(imntfSum)

# * check L4 mapped to IMN
# [1] "038 DG-PIR Ex IMN"      "045 OB-STR-CTX Inh IMN"
imn_sc_ids <- c(38, 45)
doesMappedIMN <- vapply(imn.top3sc, \(x) {
  if (any(as.integer(names(x)) %in% imn_sc_ids)) {
    return(TRUE)
  }
  return(FALSE)
}, TRUE)
L4mappedIMN <- names(imn.top3sc)[doesMappedIMN]
imn.top3sc[L4mappedIMN]
as.tibble(imn.top3sc[L4mappedIMN]) |> print(20)
t <- vapply(seq_along(L4mappedIMN), \(i) {
  l4 <- L4mappedIMN[i]
  tops <- imn.top3sc[[l4]]
  paste(l4, paste(names(tops), round(tops,2), sep = ":", collapse = ","),
    sep = "\t")
}, "1\t1:0.20")
as.tibble(t) |> print(n=20)
##  1 "21-1-2-1\t45:0.39,41:0.31"         
##  2 "12-2-3-5\t45:0.62"                 
##  3 "12-2-3-3\t39:0.67,45:0.17"         
##  4 "21-3-2-1\t45:0.75"                 
##  5 "21-3-1-1\t45:0.63,42:0.09"         
##  6 "36-2-2-1\t38:0.63,45:0.07"         
##  7 "7-5-4-1\t9:0.74,38:0.06"           
##  8 "21-2-2-1\t38:0.47,45:0.29"         
##  9 "21-2-1-1\t38:0.47,45:0.23"         
## 10 "21-6-2-1\t40:0.68,45:0.19"         
## 11 "36-1-1-1\t38:0.52,45:0.34"         
## 12 "7-11-5-1\t262:0.12,38:0.11,36:0.07"
## 13 "36-1-2-1\t38:0.33,45:0.29"         
## 14 "36-2-1-1\t38:0.77,45:0.11"         
## 15 "15-1-5-1\t38:0.45,37:0.44"         

# * prepare IMN to all the other L4 scores
imn2L4.scMat <- imn.scMat[as.integer(rownames(imn.scMat)) %in% imn_sc_ids, ]
saveRDS(imn2L4.scMat, file.path(outdir, "IMNtoL4_subclassMat.rds"))
