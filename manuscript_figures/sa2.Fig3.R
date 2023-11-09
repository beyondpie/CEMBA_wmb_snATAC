source("paper.R")
# * configs
atacMeta <- readRDS(cembav2env$sa2metaFile)
peakBed <- Sa2PeakCalling$loadFinalPeakBed()
scnmMap <- cembav2env$load.sa2.subclass.names()
nds <- 5000

# get class_color, nt_color, subclass to class and nt
sc2class <- unique(
  atacMeta[,
    c("subclass_label_v3", "class_label_v3", "class_color_v3",
      "NT_v3", "nt_type_color_v3")])
NN_type_color <- "grey"
sc2class$nt_type_color_v3[
  is.na(sc2class$nt_type_color_v3)] <- NN_type_color
rownames(scnmMap) <- scnmMap$subclass_label
sc2class$subclass <- scnmMap[sc2class$subclass_label_v3, "subclass_label_peak"]
rownames(sc2class) <- sc2class$subclass
class2color <- unique(sc2class[, c("class_label_v3", "class_color_v3")])
class_color <- class2color$class_color_v3
names(class_color) <- class2color$class_label_v3
nt2color <- unique(sc2class[, c("NT_v3", "nt_type_color_v3")])
nt_color <- nt2color$nt_type_color_v3
names(nt_color) <- nt2color$NT_v3

# get subclass in order for ATAC-seq
# * load atac signals
logCPM.scbyp <- Sa2PeakCalling$read.sa2CPM.pbysc(
  cap = 0.99, runlog1p = TRUE, scbyp = TRUE)
# scs_ord: ordered subclasses for ATAC-seq
subclassATAC <- rownames(logCPM.scbyp)
scs_id <- with(scnmMap,
  subclass_id[match(subclassATAC, subclass_label_peak)])
scs_ord <- subclassATAC[order(scs_id)]

# get subclass to number of cell
sc2n <- table(atacMeta$subclass_label_v3) |>
  as.data.frame() |>
  setNames(object = _, nm = c("subclass", "ncell")) |>
  mutate(subclass = levels(subclass)[subclass]) |>
  mutate(subclass = scnmMap[match(subclass, scnmMap$subclass_label),
    "subclass_label_peak"])
rownames(sc2n) <- sc2n$subclass

# ======Fig.3c=========
# Heatmap annotations
# * set subclass annotation
subclass_annot <- ComplexHeatmap::HeatmapAnnotation(
  class = sc2class[scs_ord, "class_label_v3"],
  col = list(class = class_color),
  which = "column",
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_label = "class",
  annotation_legend_param = list(
    class = list(
      nrow = 6,
      title_position = "leftcenter-rot")
  )
)
nt_annot <- ComplexHeatmap::HeatmapAnnotation(
  nt = sc2class[scs_ord, "NT_v3"],
  col = list(nt = nt_color),
  which = "column",
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_label = "nt",
  annotation_legend_param = list(
    nt = list(
      nrow = 2,
      title_position = "topcenter",
      title = "NT type")
  )
)
# add subclass annot text
annotIndex1 <- c(
  seq(from = 1, to = 100, by = 6),
  seq(from = 254, to = 275, by = 2))
scnmAnnot1 <- ComplexHeatmap::columnAnnotation(
  subclass = anno_mark(
    at = annotIndex1,
    labels = scs_ord[annotIndex1],
    labels_gp = grid::gpar(fontsize = 9))
)
annotIndex2 <- seq(from = 100, to = 253, by = 6)
scnmAnnot2 <- ComplexHeatmap::columnAnnotation(
  subclass = anno_mark(
    at = annotIndex2,
    labels = scs_ord[annotIndex2],
    labels_gp = grid::gpar(fontsize = 9))
)

# motif color sets
## Before filtering low motifs:  428 motifs.
## After filtering low motifs:  375 motifs.
motifColor <- circlize::colorRamp2(
  c(0, 10, 50, 100), c("white", "gray90", "gray70", "red"))
lgdMotif <- list(
  title = "-log10(p-value)",
  at = c(10, 100),
  labels = c(10, 100),
  direction = "horizontal"
)


# load scRNA-seq
RNAgbysc <- readRDS(file.path(projdir,
    "04.cCREgene/sa2.cicero/out/sa2pdcsum",
    "sa2.allen.avg.logCPM.gbysc.rds"))

# load ppdc
ppdc.with.cor <- data.table::fread(
  file = file.path(projdir, "04.cCREgene/sa2.cicero/out/sa2pdcsum",
    "mba.whole.sa2subclass.pearson.pos.pdc.alignv1.tsv"),
  header = TRUE, sep = "\t", data.table = FALSE)
genes <- vapply(ppdc.with.cor$conns, function(i) {
  unlist(strsplit(i, split = "|", fixed = TRUE))[1]
}, FUN.VALUE = "a")
distals <- vapply(ppdc.with.cor$conns, function(i) {
  unlist(strsplit(i, split = "|", fixed = TRUE))[2]
}, FUN.VALUE = "peak1")
ppdc.with.cor$gene <- genes
ppdc.with.cor$distal <- distals
distal2gene <- ppdc.with.cor |>
  group_by(distal) |>
  summarise(geneByPeak = gene[which.max(pcc)]) |>
  as.data.frame()
rownames(distal2gene) <- distal2gene$distal

# * load motif
nmod <- 54
# 428 motif in total
motifList <- lapply(seq_len(nmod), function(r) {
  a <- data.table::fread(file = file.path(projdir, "figures/data",
    "fig3/nmf.ppdcpeak.r54.motif",
    paste0("homer_n", r), "knownResults.txt"), sep = "\t",
    header = TRUE, data.table = FALSE)
  motifnms <- vapply(a[,1], function(i) {
    unlist(strsplit(i, split = "/", fixed = TRUE))[1]
  }, FUN.VALUE = "CTCF")
  logp <- a[, 4]
  motifnms <- vapply(a[,1], function(i) {
    unlist(strsplit(i, split = "/", fixed = TRUE))[1]
  }, FUN.VALUE = "CTCF")
  result <- data.frame(motif = motifnms, logp = logp)
  return(result)
})
names(motifList) <- paste0("m", seq_len(nmod))
all.motif.nms <- unique(unlist(
  lapply(motifList, function(i) {return(i$motif)})
))
logp <- (-10)
motifs.filtered <- lapply(motifList, \(m) {
  with(m, motif[logp <= logp])
})

length(unique(unlist(motifs.filtered)))

# load NMF result
ppdcStatPeak <- loadStatPeak.NMF(
  file.path(projdir, "figures/data",
    "fig3/sa2_ppdcpeak_nmf", "nmfPmat.ppdcpeak.r54.n0.statW")
)
nmod <- length(unique(ppdcStatPeak$moduleN))
statPeak.ds <- downsampleStatPeak.1(
  statPeak = ppdcStatPeak,
  mod.ordered = NULL,
  size = floor(nds / nmod),
  seed = 2022)
peak.ds <- statPeak.ds$peak

# * load atac NMF
ordModules <- orderModules.cpm.statPeak(
  cpm = logCPM.scbyp[scs_ord, statPeak.ds$peak],
  statPeak = statPeak.ds)

data.table::fwrite(
  x = data.frame(mod.order = ordModules),
  file = "sa2.ppdc.nnf.module.order.txt",
  col.names = FALSE
)
statPeak.ods <- orderStatPeakByModuleRank(
  statPeak = statPeak.ds,
  mod.ordered = ordModules)
ataclogCPM <- logCPM.scbyp[scs_ord, statPeak.ods$peak]

# set atac color bar
lowATAC <- quantile(ataclogCPM, 0.005)
highATAC <- quantile(ataclogCPM, 0.995)
ataclogCPMColor <- circlize::colorRamp2(
  seq(lowATAC, highATAC, length = 60),
  viridis::viridis(60))
lgdATAC <- list(
  title = "log(CPM + 1)",
  at = c(lowATAC, highATAC),
  labels = c(round(lowATAC, 1), round(highATAC, 1)),
  direction = "horizontal")

# align allen RNA-seq
gene.ordbydistal <- with(distal2gene,
  geneByPeak[match(statPeak.ods$peak, distal)])
RNAlogCPM.ds <- RNAgbysc[gene.ordbydistal, scs_ord] |>
  scaleByRow(mat = _)

# set rna color bar
lowRNA <- quantile(RNAlogCPM.ds, 0.05)
highRNA <- quantile(RNAlogCPM.ds, 0.95)
allenRNAColor <- circlize::colorRamp2(
  seq(lowRNA, highRNA, length = 60),
  viridis::magma(60))
lgdRNA <- list(
  title = "z-score of log(CPM+1)",
  at = c(lowRNA, highRNA),
  labels = c(round(lowRNA, 1), round(highRNA, 1)),
  direction = "horizontal")

# * heatmap
hmATAC <- ComplexHeatmap::Heatmap(
  matrix = t(ataclogCPM),
  col = ataclogCPMColor,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 6),
  show_column_names = FALSE,
  column_names_gp = grid::gpar(fontsize = 5),
  top_annotation = scnmAnnot1,
  left_annotation = NULL,
  bottom_annotation = subclass_annot,
  use_raster = TRUE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = lgdATAC,
  width = 1
)

hmRNA <- ComplexHeatmap::Heatmap(
  mat = RNAlogCPM.ds,
  col = allenRNAColor,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 6),
  show_column_names = FALSE,
  top_annotation = scnmAnnot2,
  bottom_annotation = nt_annot,
  left_annotation = NULL,
  use_raster = TRUE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = lgdRNA,
  width = 1
)
## hm.ppdc <- ComplexHeatmap::draw(hmATAC + hmRNA,
##   heatmap_legend_side = "bottom",
##   merge_legend = TRUE)

# * add motif
moduleRepeats <- table(statPeak.ods$moduleN)
motifMat <- getMotifMat(
  modules = ordModules,
  motifs = all.motif.nms,
  motifList = motifList,
  globalModules = NULL,
  lowHighEnrichedModules = 1,
  highEnrichNegLogp = 10,
  moduleRepeats = NULL
)
motifColumns <- getOrderOfMotifMat(
  motifMat = motifMat, k = 40, topk = 1, capped = 100,
  enrichBound = 90, topGlobal = 5)
motifMat <- motifMat[ , motifColumns]

expandModules <- unlist(lapply(rownames(motifMat), function(i) {
  rep(i, moduleRepeats[i])
}))
expMotifMat <- motifMat[expandModules, ]

motifColor <- circlize::colorRamp2(
  c(0, 10, 50, 100),
  c("white", "gray90", "gray90", "red"))
lgdMotif <- list(
  title = latex2exp::TeX(r"(-$\log_{10}$p)"),
  at = c(10, 100),
  labels = c(10, 100),
  direction = "horizontal")

# determine motif annot
## tmp <- data.frame(
##   index = seq_len(ncol(expMotifMat)),
##   motif = colnames(expMotifMat)
## )
## write.table(tmp, file = "tmp.moti.nms.csv", quote = FALSE,sep = ",",
##   row.names = FALSE, col.names = FALSE)
pos <- sort(unique(c(
  1:5, sample(6:139, size = 15, replace = FALSE),
  33, 140, 158, 160,
  c(326, 328, 332, 335, 338:340),
  sample(141:326, size = 10, replace = FALSE),
  sample(340:375, size = 5, replace = FALSE)
)))

motif_annot <- columnAnnotation(
  motif = anno_mark(at = pos,
    labels = colnames(expMotifMat)[pos],
  labels_gp = grid::gpar(fontsize = 9))
)

hmMotif <- ComplexHeatmap::Heatmap(
  mat = expMotifMat,
  ## mat = motifMat,
  col = motifColor,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = grid::gpar(fontsize = 3),
  top_annotation = motif_annot,
  left_annotation = NULL,
  use_raster = TRUE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = lgdMotif,
  na_col = "white",
  width = 1.25
)
# * put all the heatmap together
hm.ppdc <- ComplexHeatmap::draw(hmATAC + hmRNA + hmMotif,
  heatmap_legend_side = "bottom",
  merge_legend = TRUE)
withr::with_pdf(
  new = file.path(projdir, "figures", "out/figures",
    "sa2.Fig3.ppdc.heatmap.pdf"),
  code = {print(hm.ppdc)},
  height = 10, width = 15
)

# ===========Fig.3b==============
corRmdFile <- file.path(projdir, "04.cCREgene",
  "sa2.cicero", "out/sa2pdcsum", "sa2subclass.cor.pearson.rdm.shuf.csv")
corRealFile <- file.path(projdir, "04.cCREgene",
  "sa2.cicero", "out/sa2pdcsum", "sa2subclass.cor.pearson.real.csv")
pdcCorFile <- file.path(projdir, "04.cCREgene",
  "sa2.cicero", "out/sa2pdcsum", "sa2subclass.pdc.cor.pearson.fdr.alignv1.tsv")

plotCorDensity(
  corRealFile = corRealFile,
  corRdmFile = corRmdFile,
  pdcCorFile = pdcCorFile,
  corMethod = "PCC",
  outdir = file.path(projdir, "figures",
    "out/figures"),
  prefix = "sa2.")

# ===========Supplementary Figures==========
# * all connections for each subclass
outfigdir <- file.path(projdir,
  "figures", "out", "figure3")
if(!dir.exists(outfigdir)) {
  dir.create(outfigdir)
}
connStatAll <- data.table::fread(
  file = file.path(projdir, "04.cCREgene",
    "sa2.cicero", "out", "sa2pdcsum",
    "mba.whole.sa2subclass.merge.conns.sta.all"),
  header = FALSE,
  sep = "\t",
  data.table = FALSE
) |> setNames(
  object = _,
  nm = c("cluster", "totConns", "ddConns", "ppConns", "pdConns", "geneN", "creN")
)

rownames(connStatAll) <- connStatAll$cluster
sum(connStatAll$totConns)
## [1] 56499364
sum(connStatAll$pdConns)
## [1] 16964610
connStatAll$mainclass <- sc2class[rownames(connStatAll), "class_label_v3"]
connStatAll$ncell <- sc2n[rownames(connStatAll), "ncell"]
topk <- 10
top.all <- connStatAll[
  order(connStatAll$totConns,decreasing = T)[1:topk], ]
(
  p.conn.all <- ggplot(data = connStatAll,
  aes(x = ncell, y = totConns, color = mainclass)) +
    geom_point(size= 3, alpha = 0.8) +
    ggrepel::geom_text_repel(data = top.all,
      mapping = aes(label = cluster), size = 3)+
   ## coord_trans(x="log10", y="log10") +
   scale_y_continuous(label=scientific_10, trans = 'log10', n.breaks = 10) + 
   scale_x_continuous(trans='log10', n.breaks = 10, label=scientific_10) +
   xlab("# of nuclei in cell subclass") +
   ylab("# of conns between all the cCREs")+
   scale_color_manual(values = class_color) + 
   theme_bw() +
   theme(
     axis.text.y = element_text(colour = "black", size = 14),
     axis.text.x = element_text(colour = "black", size = 13, angle = 40, hjust=1),
     axis.title = element_text(colour = "black", size = 15),
     legend.position = "bottom",
     legend.text = element_text(colour = "black", size = 11)
   )
)
ggsave(p.conn.all, filename = file.path(outfigdir,
  "sa2.allconn.cicero.per.subclass.pdf"),
  width = 10, height=  12)

topk <- 10
top.pdc <- connStatAll[
  order(connStatAll$pdConns ,decreasing = T)[1:topk], ]

(
  p.conn.pdc <- ggplot(data = connStatAll,
  aes(x = ncell, y = pdConns, color = mainclass)) +
   geom_point(size= 3, alpha = 0.8) +
    ggrepel::geom_text_repel(data = top.pdc,
      mapping = aes(label = cluster), size = 3)+
   ## coord_trans(x="log10", y="log10") +
   scale_y_continuous(label=scientific_10, trans = 'log10', n.breaks = 10) + 
   scale_x_continuous(trans='log10', n.breaks = 10, label=scientific_10) +
   xlab("# of nuclei in cell subclass") +
   ylab("# of conns between distal-proximal cCREs")+
   scale_color_manual(values = class_color) + 
   theme_bw() +
   theme(
     axis.text.y = element_text(colour = "black", size = 14),
     axis.text.x = element_text(colour = "black", size = 13, angle = 40, hjust=1),
     axis.title = element_text(colour = "black", size = 15),
     legend.position = "bottom",
     legend.text = element_text(colour = "black", size = 11)
   )
)
ggsave(p.conn.pdc, filename = file.path(outfigdir,
  "sa2.all.conn.pdc.cicero.per.subclass.pdf"),
  width = 10, height=  12)

# * histogram
pdc.dist.all <- data.table::fread(
  file = cembav2env$sa2.all.pdc.dist,
  header = TRUE,
  data.table = FALSE
)

genes <- unique(vapply(pdc.dist.all$conns, function(i){
  unlist(strsplit(i, split = "|", fixed = TRUE))[1]
}, FUN.VALUE = "Xkr4"))

cCREs <- unique(vapply(pdc.dist.all$conns, function(i){
  unlist(strsplit(i, split = "|", fixed = TRUE))[2]
}, FUN.VALUE = "peak1"))

## length(unique(pdc.dist.all$conns))
## [1] 1997184
## length(genes)
## [1] 24414
## > length(cCREs)
## [1] 502704


mdist.all <- median(pdc.dist.all$distance)
(
  hist.dist.all <- ggplot(data = pdc.dist.all, aes(x = distance)) +
    geom_histogram(fill = "lightblue", bins = 100) +
    labs(x = "Distance between proximal-distal connections", y = "Count") +
    geom_vline(aes(xintercept = mdist.all), col = 'darkblue',
      linewidth = 1) +
    scale_y_continuous(label = scales::comma) +
    scale_x_continuous(label = scales::comma) +
    theme_bw() +
    theme(
      axis.title = element_text(colour = "black", size = 15),
      axis.text.x = element_text(colour = "black",size = 12),
      axis.text.y = element_text(colour = "black",size = 12)
    ) +
    ggtitle(paste0("median distance: ",mdist.all))
)
ggsave(hist.dist.all,
  filename = file.path(outfigdir, "sa2.hist.all.pdc.pdf"),
  height = 10, width = 10)
  
ppdc <- data.table::fread(file = cembav2env$sa2.ppdc.pearson.tsv,
  header = TRUE, sep = "\t", data.table = FALSE)
dist.ppdc <- pdc.dist.all[pdc.dist.all$conns %in% ppdc$conns, ]
mdist.ppdc <- median(dist.ppdc$distance)
(hist.dist.ppdc <- ggplot(data = dist.ppdc, aes(x = distance)) +
   geom_histogram(fill = "lightblue", bins = 100) +
   labs(x = "Distance between positive proximal-distal cCREs", y = "Count") +
   geom_vline(aes(xintercept = mdist.ppdc), col = 'darkblue',
     linewidth = 1) +
   scale_y_continuous(label = scales::comma) +
   scale_x_continuous(label = scales::comma) +
   theme_bw() +
   theme(
     axis.title = element_text(colour = "black", size = 15),
     axis.text.x = element_text(colour = "black",size = 12),
     axis.text.y = element_text(colour = "black",size = 12)
   ) +
   ggtitle(paste0("median distance: ", mdist.ppdc))

)
ggsave(hist.dist.ppdc,
  filename = file.path(outfigdir, "sa2.hist.pos.pdc.pdf"),
  height = 10, width = 10)
genes <- unique(vapply(ppdc$conns, function(i){
  unlist(strsplit(i, split = "|", fixed = TRUE))[1]
}, FUN.VALUE = "Xkr4"))

cCREs <- unique(vapply(ppdc$conns, function(i){
  unlist(strsplit(i, split = "|", fixed = TRUE))[2]
}, FUN.VALUE = "peak1"))
## + + > length(genes)
## [1] 20703
## > length(cCREs)
## [1] 281200


dist.ppdc2 <- unique(pdc.dist.all[pdc.dist.all$conns %in% ppdc$conns,
  c("conns", "distance")])
mdist.ppdc2 <- median(dist.ppdc2$distance)

gene2cCREs <- t(vapply(ppdc$conns, function(i) {
  unlist(strsplit(i, split = "|", fixed = TRUE))
}, FUN.VALUE = c("Xkr4", "peak1")))
gene2cCREs <- unique(as.data.frame(gene2cCREs))
colnames(gene2cCREs) <- c("gene", "cCRE")

cCREsPerGene <- aggregate(x = gene2cCREs$cCRE,
  by = list(gene = gene2cCREs[,1]), length)
withr::with_pdf(
  new = file.path(outfigdir, "sa2.cCREpergene.ppdc.pdf"),
  code = {
    boxplot(cCREsPerGene$x, col = "lightblue", notch = FALSE,
      xlab = "# cCREs / genes", ylab = "Count",
      outline = FALSE)
  },
  width = 4, height = 10
)
boxplot(cCREsPerGene$x, col = "lightblue", notch = FALSE,
  xlab = "# cCREs / genes", ylab = "Count",
  outline = FALSE)

  
# * negative pdc
corReal <- loadCor(cembav2env$sa2.pdc.cor.pearson.real.file)
corReal <- corReal[!is.na(corReal$cor), ]
corReal.fdr <- data.table::fread(file = cembav2env$sa2.pdc.cor.pearson.fdr.tsv,
  header = TRUE, sep = "\t",
  data.table = FALSE)

(minpcc <- min(
  with(corReal.fdr, cor[(cor >= 0) & (FDR <=0.01) ])))
(negpcc <- max(
  with(corReal.fdr, cor[(cor <= 0) & (FDR <=0.01) ])))

## length(unique(corReal.fdr[corReal.fdr$cor <= negpcc, "pdc"]))
## [1] 107413


n_pos <- sum(with(corReal.fdr, (cor >=0) & (FDR <= 0.01)))
n_neg <- sum(with(corReal.fdr, (cor < 0) & (FDR <= 0.01)))
x_text <- minpcc + 0.01
y_text <- 3
textPlus <- paste(format(n_pos, big.mark = ",", scientific = FALSE),
  "positively correlated pairs")

x2_text <- negpcc-0.01
y2_text <- 3
textNeg <- paste(format(n_neg, big.mark = ",", scientific = FALSE),
  "negatively correlated pairs")

corRdmShuf <- loadCor(cembav2env$sa2.pdc.cor.pearson.rdmshuf.file)
corRdmShuf <- corRdmShuf[!is.na(corRdmShuf$cor), ]

set.seed(2023)
corRdmShuf <- corRdmShuf[sample(nrow(corRdmShuf), nrow(corReal)),]
corsum <- rbind(corReal,corRdmShuf)
corsum$class <- c(
  rep("Real", nrow(corReal)), rep("RdmShuf", nrow(corReal)))
p.density <- ggplot() +
  geom_density(data = corRdmShuf, aes(x= cor),
    fill = "lightgray", colour = NA) +
  geom_density(data = corReal, aes(x = cor),
    colour = "black") +
  geom_vline(xintercept = minpcc, linetype = "dashed") +
  geom_vline(xintercept = negpcc, linetype = "dashed") +
  geom_area(data = subset(d, x >= minpcc), aes(x = x, y = y),
    fill = "darkred", alpha = 0.4) +
  geom_area(data = subset(d, x <= negpcc), aes(x = x, y = y),
    fill = "darkblue", alpha = 0.4) +
  geom_text(mapping = aes(label = textPlus, x = x_text + 0.2,
    y = y_text, angle = 90), size = 8) + 
  geom_text(mapping = aes(label = textNeg, x = x2_text-0.2,
    y = y2_text, angle = 90), size = 8) + 
  scale_x_continuous(limits = c(-1.05,1.05), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("PCC") +
  ylab("Density") + 
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", size =18))
p.density

ggsave(p.density,
  filename = file.path(outfigdir, "sa2.neg.pos.vs.shuf.cor.pearson.pdf"),
  width = 10, height = 10)






