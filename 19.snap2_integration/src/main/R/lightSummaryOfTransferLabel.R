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
  getMaxColScoreWithName)

# * main
args <- commandArgs(trailingOnly = TRUE)
kAnchor <- as.integer(args[1])
group <- args[2]
allenTech <- args[3]
fromdir <- file.path(here::here(), "19.snap2_integration",
  "out/tf_tscc", group)
## kAnchor <- 70
feaName <- "k8"
# allenTech <- "allen"
method <- "cca"

outdir <- file.path(here::here(), "19.snap2_integration",
  "out/tf_lightsum")
if(!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

message("set up tfSum.")
tfSum <- TransferLabelSum$new(
  kAnchor = kAnchor,
  feaName = feaName,
  dsNum = 50,
  allenTech = allenTech,
  group = group,
  method = method
)

message("load transfer label data.")
tfSum$loadtfLabel(fromdir = fromdir)

message("get transfer labels.")
refCols <- c("subclass_id", "supertype_id", "cl")
queryCol <- "L4"
query2refList <- purrr::map(refCols, \(ref) {
  tfVoteMat <- tfSum$getAvgVoteMat.queryGroup(
    tfSum, refCol = ref, queryCol = queryCol)
  L4nms <- colnames(tfVoteMat)
  maxScore <- getMaxColScoreWithName(tfVoteMat)
  data.frame(
    query = L4nms,
    cl = names(maxScore),
    clScore = maxScore,
    row.names = L4nms
  ) |> setNames(
    object = _,
    nm = c(queryCol,
      paste0(ref, "_tf"),
      paste0(ref, "_tfscore")))
})
query2Refs <- Reduce(
  function(x, y) {
    merge(x, y, by = queryCol, all = TRUE)
  }, query2refList)

message("map to subclass, supertype based on AllenAnnot and cl.")
allenAnnotMeta <- Sa2Integration$loadAllenAnnot()
query2Refs$subclass_id <- allenAnnotMeta[query2Refs$cl_tf, "subclass_id"]
query2Refs$supertype_id <- allenAnnotMeta[query2Refs$cl_tf, "supertype_id"]

message("add snapatac2 clustering size info.")
barcode2L4 <- cluSumBySa2$loadbarcode2L4()
L4size <- table(barcode2L4$L4) |> toNamedArray.1dtable()
query2Refs$L4size <- L4size[query2Refs[[queryCol]]]

outfnm <- file.path(outdir,
  paste0(tfSum$tfnmprefix(), ".lightsum.rds"))
message("prepare for light summary to: ", outfnm)
saveRDS(query2Refs, outfnm)
message("Done")

