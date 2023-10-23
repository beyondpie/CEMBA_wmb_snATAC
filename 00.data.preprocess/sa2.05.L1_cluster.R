library(data.table)
library(Matrix)
library(ggplot2)

packdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env)
import::from(.from = "peak.R", .directory = packdir,
  getTopRankClusterForPeakModule,
  rankModules)
import::from(.from = "utils.R", .directory = packdir,
  colFn)

# * load meta
atacMeta <- readRDS(cembav2env$allmetaFile)
hcL2 <- readRDS(cembav2env$L2HCFile)
L2.order <- hcL2$labels[hcL2$order]
ignoredCols <- c("LQ")
L2.order <- L2.order[!L2.order %in% ignoredCols]
 

# * load SnapATAC2 L1 clustering result
L1s <- data.table::fread(
  file = "snapatac2_pp_out/eval_l1/L1_sds.csv",
  header = TRUE,
  data.table = FALSE)
rownames(L1s) <- L1s$barcodes
L1s$leiden_all <- L1s$leiden_all + 1
L1s$leiden_top <- L1s$leiden_top + 1
tbl1 <- with(L1s, table(leiden_all, leiden_top))

# plt result
r <- tbl1
r <- r / rowSums(r)
class(r) <- "matrix"
r.melt <- reshape2::melt(data = r, as.is = TRUE)
(
  dotplot.r <- ggplot(data = r.melt,
    aes(x = leiden_all, y = leiden_top)) +
    geom_point(aes(color = value), size = 3) +
    scale_color_gradient(low = "white", high = "red") +
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
    axis.text = element_text(colour= "black", size = 10),
    axis.title = element_text(colour = "black", size = 12)
    )
)

# plot consensus matrix
barcodes <- intersect(
  rownames(atacMeta), rownames(L1s)
)
sL1s <- L1s[barcodes, ]
sL1s$L2Annot2 <- atacMeta[rownames(sL1s), "L2Annot2"]
sL1s <- sL1s[!(sL1s$L2Annot2 %in% "LQ"), ]
r <- with(sL1s, table(leiden_all, L2Annot2))
# normalize by leiden_all
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
r.melt$leiden_all <- factor(r.melt$leiden_all, levels = leiden.order)

(
  dotplot.r <- ggplot(data = r.melt,
    aes(x = leiden_all, y = L2Annot2)) +
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
    axis.title = element_text(colour = "black", size = 12)
    )
)




# plot umap with labels
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

umapAllFeat <- data.table::fread(
  file = file.path("snapatac2_pp_out/eval_l1",
    "L1_UMAP_allfeat.csv"),
  header = TRUE,
  data.table = FALSE
)
umapAllFeat$cluster <- L1s[umapAllFeat$barcodes, "leiden_all"]
set.seed(2023)
n <- 30000
r1 <- umapAllFeat[sort(sample(seq_len(nrow(umapAllFeat)),
  size = n, replace = FALSE)), ]
r1$cluster <- factor(r1$cluster)
r1 <- addLabel(r1, onCol = "cluster")
(
  pUMAP.allfeat <- ggplot(data = r1, aes(x = UMAP1, y = UMAP2)) +
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

# umap top feat
umapTopFeat <- data.table::fread(
  file = file.path("snapatac2_pp_out/eval_l1",
    "L1_UMAP_topfeat.csv"),
  header = TRUE,
  data.table = FALSE
)
umapTopFeat$cluster <- L1s[umapTopFeat$barcodes, "leiden_top"]
set.seed(2023)
n <- 30000
r2 <- umapTopFeat[sort(sample(seq_len(nrow(umapTopFeat)),
  size = n, replace = FALSE)), ]
r2$cluster <- factor(r2$cluster)
r2 <- addLabel(r2, onCol = "cluster")
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

# * Annot L1s top features by L2
barcodes <- intersect(rownames(atacMeta),
  rownames(L1s))
r3 <- data.frame(
  barcode = barcodes,
  L1Top = L1s[barcodes, "leiden_top"],
  L2Annot2 = atacMeta[barcodes, "L2Annot2"]
)

r <- with(r3, table(L2Annot2, L1Top))
r <- r / rowSums(r)
class(r) <- "matrix"
r.melt <- reshape2::melt(data = r, as.is = TRUE)
(
  dotplot.r <- ggplot(data = r.melt,
    aes(x = L1Top, y = L2Annot2)) +
    geom_point(aes(color = value), size = 3) +
    scale_color_gradient(low = "white", high = "red") +
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
    axis.text = element_text(colour= "black", size = 10),
    axis.title = element_text(colour = "black", size = 12)
    )
)
rownames(umapTopFeat) <- umapTopFeat$barcodes
umap3 <- umapTopFeat[r3$barcode, ]
rownames(r3) <- r3$barcode
umap3$L2Annot2 <- r3[, "L2Annot2"]
umap3 <- umap3[sort(sample(seq_len(nrow(umap3)), n, replace = FALSE)) ,]
umap3$cluster <- factor(umap3$cluster)
umap3 <- addLabel(umap3, "L2Annot2")

(
  pUMAP.topfeat <- ggplot(data = umap3, aes(x = UMAP1, y = UMAP2)) +
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

# * Annot L1s all features by L2
barcodes <- intersect(rownames(atacMeta),
  rownames(L1s))
r3 <- data.frame(
  barcode = barcodes,
  L1Top = L1s[barcodes, "leiden_all"],
  L2Annot2 = atacMeta[barcodes, "L2Annot2"]
)

r <- with(r3, table(L2Annot2, L1Top))
r <- r / rowSums(r)
class(r) <- "matrix"
r.melt <- reshape2::melt(data = r, as.is = TRUE)

(
  dotplot.r <- ggplot(data = r.melt,
    aes(x = L1Top, y = L2Annot2)) +
    geom_point(aes(color = value), size = 3) +
    scale_color_gradient(low = "white", high = "red") +
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
    axis.text = element_text(colour= "black", size = 10),
    axis.title = element_text(colour = "black", size = 12)
    )
)

rownames(umapAllFeat) <- umapAllFeat$barcodes
umap3 <- umapAllFeat[r3$barcode, ]
rownames(r3) <- r3$barcode
umap3$L2Annot2 <- r3[, "L2Annot2"]
umap3 <- umap3[sort(sample(seq_len(nrow(umap3)), n, replace = FALSE)) ,]
umap3$cluster <- factor(umap3$cluster)
umap3 <- addLabel(umap3, "L2Annot2")

(
  pUMAP.topfeat <- ggplot(data = umap3, aes(x = UMAP1, y = UMAP2)) +
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




