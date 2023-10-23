library(data.table)
library(Matrix)
library(ggplot2)


packdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = packdir, cembav2env)
import::from(.from = "peak.R", .directory = packdir,
  getTopRankClusterForPeakModule,
  rankModules)
import::from(.from = "utils.R", .directory = packdir,
  colFn)

atacMeta <- readRDS(cembav2env$allmetaFile)
hcL2 <- readRDS(cembav2env$L2HCFile)
L2.order <- hcL2$labels[hcL2$order]
ignoredCols <- c("LQ")
L2.order <- L2.order[!L2.order %in% ignoredCols]

# main
L1s <- data.table::fread(
  file = "../result/clustering_sum_L1/sa2_L1_r0.4_barcodes2id.csv",
  header = TRUE,
  data.table= FALSE,
  sep = ","
)
rownames(L1s) <- L1s$barcode
L1s$leiden <- L1s$L1 + 1

# plot consensus matrix

barcodes <- intersect(
  rownames(atacMeta), rownames(L1s)
) # 2,204,291

sL1s <- L1s[barcodes, ]
sL1s$L2Annot2 <- atacMeta[rownames(sL1s), "L2Annot2"]
sL1s <- sL1s[!(sL1s$L2Annot2 %in% "LQ"), ]
r <- with(sL1s, table(leiden, L2Annot2))
# normalize by leiden
r <- r / rowSums(r)
# normalize by L2Annot2
## r <- t(t(r) / rowSums(t(r)))
r <- r[, L2.order]

topRank.leiden <- getTopRankClusterForPeakModule(
  cluster2moduleScore = t(r),
  modules = rownames(r),
  topk = 1
)
leiden.order <- rankModules(
  topRank.leiden,
  modules = rownames(r),
  avg.fn = colFn(min)
)

r.melt <- reshape2::melt(data = r, as.is = TRUE)
r.melt$L2Annot2 <- factor(r.melt$L2Annot2, levels = L2.order)
r.melt$leiden <- factor(r.melt$leiden, levels = leiden.order)

(
  dotplot.r <- ggplot(data = r.melt,
    aes(x = leiden, y = L2Annot2)) +
    geom_point(aes(color = value), size = 3) +
    scale_color_gradient(low = "white", high = "red") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 7),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = "NA"),
    panel.spacing = unit(0.002, "lines"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.text = element_text(colour= "black", size = 10),
    axis.title = element_text(colour = "black", size = 12),
    plot.title = element_text(colour = "black", size = 14, hjust = 0.5)
  ) +
    ggtitle("Consensus Matrix between SnapATAC2 L1 and SnapATAC L2 clusters (2.2 million cells).\n Normalized by #cell in SnapATAC2 clusters.") +
    xlab("37 L1 clusters by SnapATAC2") +
    ylab("73 L2 clusters by SnapATAC")
)

# consensus matrix normalized by #cell in SnapATAC L2 clusters
r <- with(sL1s, table(leiden, L2Annot2))
# normalize by leiden
r <- t(t(r) / rowSums(t(r)))
# normalize by L2Annot2
## r <- t(t(r) / rowSums(t(r)))
r <- r[, L2.order]

topRank.leiden <- getTopRankClusterForPeakModule(
  cluster2moduleScore = t(r),
  modules = rownames(r),
  topk = 1
)
leiden.order <- rankModules(
  topRank.leiden,
  modules = rownames(r),
  avg.fn = colFn(min)
)

r.melt <- reshape2::melt(data = r, as.is = TRUE)
r.melt$L2Annot2 <- factor(r.melt$L2Annot2, levels = L2.order)
r.melt$leiden <- factor(r.melt$leiden, levels = leiden.order)

(
  dotplot.r <- ggplot(data = r.melt,
    aes(x = leiden, y = L2Annot2)) +
    geom_point(aes(color = value), size = 3) +
    scale_color_gradient(low = "white", high = "red") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 7),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = "NA"),
    panel.spacing = unit(0.002, "lines"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.text = element_text(colour= "black", size = 10),
    axis.title = element_text(colour = "black", size = 12),
    plot.title = element_text(colour = "black", size = 14, hjust = 0.5)
  ) +
    ggtitle("Consensus Matrix between SnapATAC2 L1 and SnapATAC L2 clusters (2.2 million cells).\n Normalized by #cell in 73 SnapATAC clusters.") +
    xlab("37 L1 clusters by SnapATAC2") +
    ylab("73 L2 clusters by SnapATAC")
)

# draw UMAP and then annotate SnapATAC2 L1
addLabel <- function(df,
                     onCol = "L2Annot2") {
  uLabels <- unique(df[,onCol])
  labdfList <- lapply(uLabels, function(i) {
    index <- which(df[ ,onCol] == i)
    leftNA <- rep(NA, length(index)-1)
    data.frame(rownm = rownames(df)[index],
      x = c(median(df[index, "UMAP1"]), leftNA),
      y = c(median(df[index, "UMAP2"]), leftNA),
      label = c(i, leftNA))
  })
  r <- do.call(rbind, labdfList)
  rownames(r) <- r$rownm
  df$xlab <- r[rownames(df), "x"]
  df$ylab <- r[rownames(df), "y"]
  df$label <- r[rownames(df), "label"]
  return(df)
}

L1umap <- data.table::fread(
  file = file.path("../result/clustering_sum_L1/L1_UMAP.csv"),
  header = TRUE,
  data.table = FALSE)
rownames(L1umap) <- L1umap$barcode
L1umap$cluster <- L1s[L1umap$barcode, "leiden"]
L1umap <- L1umap[barcodes, ]
L1umap$annot <- atacMeta[L1umap$barcode, "L2Annot2"]

set.seed(2023)
n <- 30000
r2 <- L1umap[sort(sample(seq_len(nrow(L1umap)),
  size = n, replace = FALSE)), ]
r2$cluster <- factor(r2$cluster)
r2 <- addLabel(r2, onCol = "annot")
(
  pUMAP.topfeat <- ggplot(data = r2, aes(x = UMAP1, y = UMAP2)) +
    geom_point(size = 0.1, aes(color = cluster), alpha = 0.5) +
    ggrepel::geom_text_repel(aes(xlab, ylab, label = label),
      size = 5, show.legend = FALSE, na.rm = TRUE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 7),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = "NA"),
    panel.spacing = unit(0.002, "lines"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
    )
)


