library(Matrix)
library(data.table)
library(R6)
library(purrr)
library(ggplot2)
proj_dir <- here::here()
import::from(
  .from = "rknn.R",
  .directory = file.path(here::here(), "package/R"),
  SingleCellGroup,
  GroupAnchor,
  initSingleCellGroup,
  getNeighborMat,
  getReciprocalPairs,
  initGroupAnchor
  )


setMajorRegion <- function(scGroup, regionMeta) {}
# orderGroup <- function(simdf, rowOrder){}
filterGroupAnchorByRegion <- function(groupAnchors){}

plotSimMatWithAnchors <- function(groupAnchors, rightOrder,
                                  dotsize = 3) {
  r <- groupAnchors$pair
  r$right <- factor(r$right, levels = rightOrder)
  r <- r[order(r$right), ]
  r$left <- factor(r$left, levels = unique(r$left))
  p <- ggplot(data = r) +
    aes(x = left, y = right) +
    geom_point(aes(color = score), size = dotsize)
  return(p)
}

plotSimMat <- function(simdf, groupAnchors,
                       rightOrder,
                       dotsize = 3){
  r <- groupAnchors$pair
  r$right <- factor(r$right, levels = rightOrder)
  r <- r[order(r$right), ]
  leftOrder <- unique(r$left)
  rr <- reshape2::melt(data = simdf, as.is =TRUE)
  colnames(rr) <- c("left", "right", "score")
  rr$left <- factor(rr$left, levels = leftOrder)
  rr$right <- factor(rr$right, levels = rightOrder)
  p <- ggplot(data = rr) +
    aes(x = left, y = right) +
    geom_point(aes(color = score), size = dotsize)
  return(p)
}

# * configs
projdir <- here::here()
nnAllenFile <- file.path(projdir, "19.snap2_integration",
  "src/main/resource",
  "nn_male_10xv3", "nn_male_10xv3_tf_by_scid_avg_lognorm.csv")
neuronAllenFile <- file.path(projdir, "19.snap2_integration",
  "src/main/resource", "neuron_male_10xv3",
  "neuron_male_10xv3_tf_by_scid_avg_lognorm.csv")
snapFile <- file.path(projdir, "19.snap2_integration",
  "src/main/resource", "sa2_gmat",
  "snap_all_tf_by_scid_avg_lognorm.csv")

# * main
theme <- theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 7),
    ## legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = "NA"),
    panel.spacing = unit(0.002, "lines"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.text.x = element_text(colour = "black", angle = 90, size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.title = element_text(colour = "black", size = 14),
    plot.title = element_text(colour = "black", face = "bold", size = 15,
      hjust = 0.5)
)

## load allen's data
nnAllen <- readPandasCSV(nnAllenFile)
neuronAllen <- readPandasCSV(neuronAllenFile) |>  _[rownames(nnAllen), ]
tfAllen <- cbind(nnAllen, neuronAllen)
orderOfId <- order(as.integer(colnames(tfAllen)))
tfAllen <- tfAllen[, orderOfId]

## load snap's data
tfSnap <- readPandasCSV(snapFile) |> _[rownames(nnAllen), ]

# * calculate distances / closeness
corMethod <- "spearman"
corSnap2Allen <- cor(tfSnap, tfAllen, method = corMethod)
groupAnchorTFSpearman <- initGroupAnchor(
  corSnap2Allen, topk = 3,
  pairLevel = c("L2", "subclass_id"),
  method = "spearman",
  featureUsed = "TF"
)

(
  pTFSpearmanWithAnchor <- plotSimMatWithAnchors(groupAnchorTFSpearman,
    colnames(corSnap2Allen), 1) +
    scale_color_gradient(low = "white", high = "red",
      limits = c(0, 0.9), na.value = "white"
    ) +
    theme +
    xlab("SnapATAC2 248 L2-level clusters") +
    ylab("Allen scRNA-seq 10Xv3 334 subclass-level clusters") +
    ggtitle("Reciprocal KNN (K=3): SnapATAC L2 - Allen Subclass")
)
ggsave(plot = pTFSpearmanWithAnchor, filename = file.path(
  projdir, "19.snap2_integration", "out/figures",
  "L2_Subclass_TF_Spearman_Anchors.pdf"
))  
  
(
  pTFSpearman <- plotSimMat(corSnap2Allen,
    groupAnchorTFSpearman,
    colnames(corSnap2Allen), 1) +
    scale_color_gradient(low = "white", high = "red",
      limits = c(0.7, 0.9), na.value = "white"
    ) +
    theme +
    xlab("SnapATAC2 248 L2-level clusters") +
    ylab("Allen scRNA-seq 10Xv3 334 subclass-level clusters") +
    ggtitle("Spearman corr: SnapATAC L2 - Allen Subclass")
)
ggsave(plot = pTFSpearman, filename = file.path(
  projdir, "19.snap2_integration", "out/figures",
  "L2_Subclass_TF_Spearman.pdf"
))

# TODO: set different K to check the mapped subclass
# TODO: use different features to check the mapped subclass
# TODO: use different Allen's data to check the mapped subclass
# TODO: add region-check

# * test
a <- initSingleCellGroup(corSnap2Allen, groupId = "1-21",
  topk = 3)
b <- initSingleCellGroup(t(corSnap2Allen), groupId = "327",
  topk = 3)
simdf <- corSnap2Allen
topk <- 3
rowscGroups <- map(rownames(simdf), initSingleCellGroup,
  simdf = simdf, topk = topk)
names(rowscGroups) <- rownames(simdf)

colscGroups <- map(colnames(simdf), initSingleCellGroup,
  simdf = t(simdf), topk = topk)
names(colscGroups) <- colnames(simdf)

pairdf <- getReciprocalPairs(rowscGroups, colscGroups, corSnap2Allen)

