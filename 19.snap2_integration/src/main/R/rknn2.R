library(data.table)
library(purrr)
library(ggplot2)
library(stringr)
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(
  .from = "rknn.R", .directory = rdir,
  SingleCellGroup,
  GroupAnchor,
  initSingleCellGroup,
  getNeighborMat,
  getReciprocalPairs,
  sumGroupAnchor,
  plotSimMatLimitedByAnchorSet)
import::from(
  .from = "utils.R", .directory = rdir,
  readPandasCSV)
import::from(
  .from = "cembav2env.R", .directory = rdir,
  Sa2Integration)
import::from(
  .from = "ggtheme.R",
  .directory = rdir,
  dotplotTheme)
import::from(
    .from = "utils.R",
    .directory = file.path(here::here(), "package/R"),
    getIndexOfUniqueName)

# * functions
sumGroupAnchor <- function(ga) {
  # how many pairs
  np <- nrow(ga$pair)
  # how many sa2
  nSa2 <- length(unique(ga$pair[, 1]))
  # how many Allen
  nAllen <- length(unique(ga$pair[, 2]))
  # how many 1-to-1
  indexUniqSa2 <- getIndexOfUniqueName(ga$pair[, 1])
  indexUniqAllen <- getIndexOfUniqueName(ga$pair[, 2])
  nUniqPair <- length(intersect(indexUniqSa2, indexUniqAllen))
  data.frame(
    key = ga$id,
    nPair = np,
    nSa2 = nSa2,
    nUniqSa2 = length(indexUniqSa2),
    nAllen = nAllen,
    nUniqAllen = length(indexUniqAllen),
    nUniqPair = nUniqPair
  )
}


# * merge allen's data
#' @return data.frame, gene by group
attachAllenAvgExp <- function(tech = "10xv3",
                              group = "subclass_id",
                              fromdir = "../resource/rknn",
                              ifSaveAsRDS = TRUE) {
  nnfnm <- file.path(fromdir,
    str_glue("nn_{tech}_{group}.avgexp.csv"))
  neuronfnm <- file.path(fromdir,
    str_glue("neuron_{tech}_{group}.avgexp.csv"))
  nndf <- readPandasCSV(nnfnm)
  neurondf <- readPandasCSV(neuronfnm) |> _[rownames(nndf), ]
  mergedf <- cbind(nndf, neurondf)
  orderOfId <- order(as.integer(colnames(mergedf)))
  r <- mergedf[, orderOfId]
  if (ifSaveAsRDS) {
    message("save data.frame into rds format.")
    saveRDS(object = r,
      file = file.path(fromdir,
        str_glue("all_{tech}_{group}.allen.avgexp.rds")))
  }
  return(r)
}

attachAllenAvgExp2 <- function(group = "subclass_id",
                               fromdir = "../resource/rknn",
                               ifSaveAsRDS = TRUE) {
  mergedf <- readPandasCSV(file.path(fromdir,
    str_glue("all_10x_multiome_{group}.avgexp.csv")))
  orderOfId <- order(as.integer(colnames(mergedf)))
  r <- mergedf[, orderOfId]
  if (ifSaveAsRDS) {
    message("save data.frame into rds format.")
    saveRDS(object = r,
      file = file.path(fromdir,
        str_glue("all_{group}.allen.avgexp.rds")))
  }
  return(r)
}


# 1. 10X v3
allen10xv3Subclass <- attachAllenAvgExp(
  tech = "10xv3", group = "subclass_id", ifSaveAsRDS = FALSE)
allen10xv3Class <- attachAllenAvgExp(
  tech = "10xv3", group = "class_id", ifSaveAsRDS = TRUE)
allen10xv3Supertype <- attachAllenAvgExp(
  tech = "10xv3", group = "supertype_id", ifSaveAsRDS = TRUE)
allen10xv3Cluster <- attachAllenAvgExp(
  tech = "10xv3", group = "cluster_id", ifSaveAsRDS = TRUE)
# 2. 10X v2
## allen10xv2 <- attachAllenAvgExp(
##   tech = "10xv2", group = "subclass_id"
## )
# 10Xv2 + 10Xv3 + multiome
allenAllSubclass <- attachAllenAvgExp2(group = "subclass_id", ifSaveAsRDS = TRUE)
allenAllClass <- attachAllenAvgExp2(group = "class_id", ifSaveAsRDS = TRUE)
allenAllSupertype <- attachAllenAvgExp2(group = "supertype_id", ifSaveAsRDS = TRUE)
allenAllCluster <- attachAllenAvgExp2(group = "cluster_id", ifSaveAsRDS = TRUE)



# * load snap data
avgsnapL2File <- file.path("../resource",
  "sa2_gmat", "sa2_sa2default_L2.avgexp.csv")
avgsnapL2 <- readPandasCSV(avgsnapL2File)

avgsnapL3File <- file.path("../resource",
  "sa2_gmat", "sa2_sa2default_L3.avgexp.csv")
avgsnapL3 <- readPandasCSV(avgsnapL3File)
# need to filter genes since not all the genes
# are in snap.



# * set up parameters
## allenAvgList <- list(
##   "class" = allen10xv3Class,
##   "subclass" = allen10xv3Subclass,
##   "supertype" = allen10xv3Supertype,
##   "cluster" = allen10xv3Cluster
## )
allenAvgList <- list(
  "class" = allenAllClass,
  "subclass" = allenAllSubclass,
  "supertype" = allenAllSupertype,
  "cluster" = allenAllCluster
  )

sa2AvgList <- list("L2" = avgsnapL2, "L3"= avgsnapL3)
geneList <- Sa2Integration$getMarkersList()
ks <- c(3, 5, 8, 10, 15)
feaList <- c("TF", "k8", "markers", "merfishmarkers")
cormethods <- c("pearson", "spearman")
mytheme <- dotplotTheme(legend.pos = "right")
outDir <- file.path(here::here(),
  "19.snap2_integration", "out", "10x_multiome")


# * rrknn
rrknn <- function(fea, corMethod, topk,
                  sa2Level = "L2",
                  allenLevel = "subclass") {
  message("rrknn: ",
    paste("allen:", allenLevel), paste("sa2:", sa2Level))
  message("Parameters: ", paste("k:", topk),
    paste("fea:", fea), paste("cor:", corMethod))
  allenMat <- allenAvgList[[allenLevel]]
  snapMat <- sa2AvgList[[sa2Level]]
  rawfeas <- geneList[[fea]]
  fInAllen <- intersect(rawfeas, rownames(allenMat))
  fInSnap <- intersect(rawfeas, rownames(snapMat))
  message("Raw features: ", length(rawfeas))
  message("Found in Allen: ", length(fInAllen))
  message("Found in Snap2: ", length(fInSnap))
  finalfeas <- intersect(fInAllen, fInSnap)
  message("final features: ", length(finalfeas))
  message("Calculate similiary/correlation: ", corMethod)
  corSnap2Allen <- cor(snapMat[finalfeas, ], allenMat[finalfeas, ],
    method = corMethod)
  message("get anchors by top ", topk)
  ga <- GroupAnchor$new(
    corSnap2Allen, topk = topk, pairLevel = c(sa2Level, allenLevel),
    method = corMethod, featureUsed = fea
  )
  saveRDS(ga, file = file.path(outDir,
    "rknn", str_glue("{ga$id}.rds")))
  message("summarize GroupAnchor.")
  sumga <- sumGroupAnchor(ga)
  print(sumga)
  saveRDS(sumga, file = file.path(outDir,
    "rknn", str_glue("{ga$id}_sum.rds")))
  sim.p <- plotSimMatLimitedByAnchorSet(ga = ga, dotsize = 1)
  lowSimScore <- quantile(ga$simSa2Allen, 0.8)
  highSimScore <- max(ga$simSa2Allen)
  nSa2Used <- ga$pair$sa2 |> sort() |> unique() |> length()
  nAllenUsed <- ga$pair$allen |> sort() |> unique() |> length()
  sim.p <- sim.p +
    mytheme +
    scale_color_gradient(low = "white", high = "red",
      limits = c(lowSimScore, highSimScore), na.value = "white") +
    xlab(str_glue("SnapATAC2 {sa2Level}-level {nSa2Used} clusters.")) +
    ylab(str_glue("Allen {allenLevel}-level {nAllenUsed} clusters.")) +
    ggtitle(ga$id)
  ggsave(plot = sim.p,
    filename = file.path(outDir, "figures",
      str_glue("{sa2Level}-{allenLevel}_{fea}_{corMethod}.pdf")),
    width = 10, height = 10)
  return(sumga)
}

rrknn2 <- function(sa2Level = "L2", allenLevel = "subclass") {
  tmpGA <- rrknn("TF", "spearman", 3)
  sumGAList <- rep(list(tmpGA), length(feaList) * length(cormethods) * length(ks))
  i <- 1
  for (fea in feaList) {
    for (m in cormethods) {
      for (k in ks) {
        ga <- rrknn(fea = fea, corMethod = m, topk = k,
          sa2Level = sa2Level, allenLevel = allenLevel)
        sumGAList[[i]] <- ga
        i <- i + 1
      }
    }
  }
  sumGA <- do.call(rbind, sumGAList)
  ## rownames(sumGA) <- sumGA$key
  write.table(sumGA,
    file = file.path(outDir, "rknn",
      str_glue("{sa2Level}-{allenLevel}.sum.rknn.csv")),
    quote = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)
  return(sumGA)
}

sumGAL2 <- rrknn2(sa2Level = "L2")
sumGAL3 <- rrknn2(sa2Level = "L3")


sumGAllenL3 <- map(c("class", "subclass", "supertype", "cluster"), \(nm) {
  rrknn2(sa2Level = "L3", allenLevel = nm)
})
names(sumGAllenL3) <- c("class", "supertype", "cluster")
saveRDS(sumGAllenL3, file.path(outDir, "rknn", "allenLevels2Sa2L3.rds"))

sumGAllenL2 <- map(c("class", "subclass", "supertype", "cluster"), \(nm) {
  rrknn2(sa2Level = "L3", allenLevel = nm)
})
names(sumGAllenL2) <- c("class", "supertype", "cluster")
saveRDS(sumGAllenL2, file.path(outDir, "rknn", "allenLevels2Sa2L2.rds"))

g <- "supertype"
sumGAllenL2[[g]][which.max(sumGAllenL2[[g]]$nAllen), ]
g <- "cluster"
sumGAllenL2[[g]][which.max(sumGAllenL2[[g]]$nAllen), ]

# * test
## fea <- "TF"
## k <- 3
## m <- "spearman"
## rrknn(fea, m, k)
