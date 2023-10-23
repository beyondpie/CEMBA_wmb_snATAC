library(data.table)
library(Matrix)
library(ggplot2)
library(stringr)
library(ggpubr)

packdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env, cluSumBySa2,
  loadL2Resos)
import::from(.from = "sa2Leiden.R", .directory = packdir,
  LeidenSum,
  newLeidenSumFromPickleFile,
  getCluster2Size,
  getBarcode2ClusterName)
import::from(.from = "utils.R", .directory = packdir,
  colFn)
# TODO: merge these common functions to utils
import::from(.from = "peak.R", .directory = packdir,
  getTopRankClusterForPeakModule,
  rankModules)
import::from(.from = "colors.R", .directory = packdir,
  SnapATACPalette)

# * load sa2 L2 clustering data
sa2L2SumDir <- file.path(here::here(),
  "17.snapatac2/result", "L2_sa2_dlt2_sum")
nm_prefix <- "sa2_clustering"
cll <- "L1"
leidenSums <- lapply(seq(0,36), function(i) {
  message("Loading cluster:", i, "'s LeisenSum.")
  f <- file.path(sa2L2SumDir, str_glue("{nm_prefix}_{cll}_{i}.pkl"))
  newLeidenSumFromPickleFile(pickleFile = f,
    cll = cll, cid = i, nm_prefix = nm_prefix,
    cidIsInt = TRUE,
    cidFromZero = TRUE,
    python = "~/mambaforge/bin/python")
})
names(leidenSums) <- vapply(leidenSums, function(i){i$cid},"0")

# * summarize the resolution results
# done in cluSumBySa2 in cemba2env.R
resos4L2 <- loadL2Reos()

# check data
getnclu <- function(leidenSum, r= 0.1) {
  j <- which(leidenSum$rs == r)
  length(table(leidenSum$label[, j]))
}

ncluFromSum <- vapply(names(leidenSums), function(nm) {
  s <- leidenSums[[nm]]
  getnclu(s, r = resos4L2[nm, "reso"])
}, FUN.VALUE = 0)

all(resos4L2$nclu == ncluFromSum)

# * load sa1 clustering result
atacMeta <- readRDS(cembav2env$allmetaFile)
hcL2 <- readRDS(cembav2env$L2HCFile)
L2.order <- hcL2$labels[hcL2$order]
ignoredCols <- c("LQ")
L2.order <- L2.order[!L2.order %in% ignoredCols]
barcodesL1 <- rownames(atacMeta)
## L3OrderFile = file.path(here::here(), "meta",
##   "L3toAllen.hc.order.txt"),

# * draw consensus matrix
barcode2Clu <- function(leidenSum, r = 0.1, i = 1) {
  ith <- which(leidenSum$rs == r)
  data.frame(
    barcode = leidenSum$barcodes,
    clu = paste(i, leidenSum$label[ , ith], sep = "_")
  )
}
barcode2CluL2s <- lapply(seq_along(leidenSums), function(i){
  s <- leidenSums[[i]]
  barcode2Clu(s, r = resos4L2[i, "reso"], i = i)
})
barcode2CluL2 <- do.call(rbind, barcode2CluL2s)
rownames(barcode2CluL2) <- barcode2CluL2$barcode
# ** L2 level
barcodes <- intersect(
  rownames(atacMeta), rownames(barcode2CluL2)
) # 2,204,291

sL1s <- barcode2CluL2[barcodes, ]
sL1s$L2Annot2 <- atacMeta[barcodes, "L2Annot2"]
sL1s <- sL1s[!(sL1s$L2Annot2 %in% "LQ"), ]

r <- with(sL1s, table(clu, L2Annot2))
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
r.melt$leiden <- factor(r.melt$clu, levels = )

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

# ** allen and L3 
Int.L3.meta <- read.table(
  file = cembav2env$IntL3MetaFile,
  sep = ",",
  header = TRUE
)
rownames(Int.L3.meta) <- Int.L3.meta$L3
Int.L3.meta$Id <- seq_len(nrow(Int.L3.meta))

a <- readRDS("../../figures/L3.hc.meta.rds")
hc.sum <- a[[4]]
sc.order.hc.class <- hc.sum$hc$labels[hc.sum$hc$order]
allen.order <- Int.L3.meta[sc.order.hc.class, "L3AllenAnnot"]

sL1s <- barcode2CluL2[barcodes, ]
sL1s$L2Annot2 <- atacMeta[barcodes, "L3"]
sL1s <- sL1s[!(sL1s$L2Annot2 %in% "LQ"), ]

r <- with(sL1s, table(clu, L2Annot2))
# normalize by leiden
## r <- r / rowSums(r)
# normalize by L2Annot2
r <- t(t(r) / rowSums(t(r)))
r <- r[, sc.order.hc.class]

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
r.melt$L2Annot2 <- factor(r.melt$L2Annot2, levels = sc.order.hc.class)
r.melt$leiden <- factor(r.melt$clu, levels = )

(
  dotplot.r <- ggplot(data = r.melt,
    aes(x = leiden, y = L2Annot2)) +
    geom_point(aes(color = value), size = 1) +
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
    axis.text = element_text(colour= "black", size = 6),
    axis.title = element_text(colour = "black", size = 6),
    plot.title = element_text(colour = "black", size = 14, hjust = 0.5)
  ) +
    ggtitle("Consensus Matrix between SnapATAC2 L1 and SnapATAC L2 clusters (2.2 million cells).\n Normalized by #cell in SnapATAC2 clusters.") +
    xlab("37 L1 clusters by SnapATAC2") +
    ylab("73 L2 clusters by SnapATAC")
)

# * prepare L3 clustering
## check barcodes order
drawUMAP <- function(leidenSum, r = 0.1, n = 20000, id = 0) {
  xy <- leidenSum$umap
  ith <- which(a$rs == r)
  l <- leidenSum$label[, ith]
  nclu <- length(unique(l))
  ncell <- nrow(xy)
  message(str_glue("{ncell} cells ",
    "{nclu} clusters under r: {r}"))
  d <- data.frame(
    x = xy[, 1],
    y = xy[, 2],
    label = factor(l, levels = seq(1, length(unique(l))))
  )
  if (nrow(d) > n){
    rows <- sample(seq_len(nrow(d)), size = n, replace = FALSE)
    d <- d[rows, ]
  }
  ggplot(data = d, mapping = aes(x = x, y = y)) +
    geom_point(size = 0.8, aes(color = label), alpha = 0.5) +
    scale_color_manual(values = SnapATACPalette) + 
    theme_minimal() +
    ggtitle(str_glue("L1-{id}: n-{ncell},r-{r},#c{nclu}")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      axis.title = element_blank()
    )
}

gplots <- lapply(seq_along(leidenSums), function(i) {
  l <- leidenSums[[i]]
  r <- resos4L2[i, "reso"]
  drawUMAP(l, r, id = i)
})


plot_all <- ggarrange(plotlist = gplots,
  ncol = 6, nrow = 7)
ggsave(filename = "../result/sa2_dlt2_L2.pdf",
  plot = plot_all, width = 15, height = 17)

## prepare L3 clustering meta data
# cluster index as 1-2: L1 Cluster Id 1's sub-cluster Id 2
# index start from 1
# - cluster2size
cluster2size <- getCluster2Size(
  leidenSumList = leidenSums,
  reso4Clu = resos4L2
)
# - barcode2id
barcode2id <- getBarcode2ClusterName(
  leidenSumList = leidenSums,
  reso4Clu = resos4L2,
  cll = "L2"
)
saveRDS(object =  list(cluster2size, barcode2id),
  file = "prepae.L3.clustering.rds")

# check the data we generated
a <- readRDS("prepae.L3.clustering.rds")
cluster2size <- a[[1]]
barcode2id <- a[[2]]

length(unique(cluster2size$clu)) == nrow(cluster2size)
rownames(cluster2size) <- cluster2size$clu
sum(cluster2size$count) == nrow(barcode2id)
length(unique(barcode2id$barcode)) == nrow(barcode2id)

t <- table(barcode2id$L2)[rownames(cluster2size)]
all(t == cluster2size$count)

all(vapply(seq_len(nrow(resos4L2)), function(i){
  cluster <- resos4L2[i, "cluster"]
  n1 <- resos4L2[i, "nclu"]
  n2 <- sum(grepl(str_glue("^{cluster}-"), cluster2size$clu))
  message(str_glue("check cluster: {cluster}, n1 {n1} vs n2 {n2}."))
  n1 == n2
}, TRUE))

# * now save files to txt and prepare L3.
write.table(x = barcode2id, file = "../resource/sa2_dlt2_L2_barcode2id.csv",
  quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
write.table(x = cluster2size, file = "../resource/sa2_dlt2_L2_cluster2size.csv",
  quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
