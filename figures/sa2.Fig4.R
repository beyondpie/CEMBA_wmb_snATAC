source("paper.R")


# * configs
outdir <- file.path(projdir, "figures", "data/fig4")
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

TF.df.co <- data.table::fread(
  file = file.path(projdir, "22.sa2GRN", "src/main/resource",
    "CisBP_ver2_Mus_musculus.motif2factors.txt"),
  header = TRUE,
  data.table = FALSE
)

# * load GRN to sparse matrix 3-column format
GRNDir <- file.path(projdir, "22.sa2GRN",
  "out", "GRN", "baseGRN_subclass")
GRNFiles <- list.files(path = GRNDir,
  all.files = FALSE, pattern = "^rawGRN")
subclasses <- vapply(GRNFiles, function(i) {
  unlist(strsplit(i, split = ".", fixed = TRUE))[2]
}, FUN.VALUE = "ZI_Pax6_Gaba")

# 267 GRNs
rawGRNs <- lapply(GRNFiles,
  loadRawGRNFile, indir = GRNDir)
names(rawGRNs) <- subclasses
## saveRDS(rawGRNs,
  ## file.path(outdir, "sa2.rawGRNs.baseGRN_subclass.rds"))

# * CellOracle default filtering
# here nlogp is -log10(p-vlaue)
lownlogp <- 3
topk <- 10000
GRNs <- lapply(rawGRNs, function(i) {
  absCut <- sort(abs(i$coef), decreasing = TRUE)[
    min(nrow(i), topk)]
  subset(i, (nlogp >= lownlogp) & (abs(coef) >= absCut))
})
names(GRNs) <- names(rawGRNs)
## saveRDS(GRNs,
##   file.path(outdir, "sa2.GRNs.baseGRN_subclass.nlogp3_topk10000.rds"))
GRNs <- readRDS(file.path(outdir,
  "sa2.GRNs.baseGRN_subclass.nlogp3_topk10000.rds"))
# 401
allTF.source <- unique(unlist(lapply(GRNs, \(i) { i$source })))
# 2985
allGene.target <- unique(unlist(lapply(GRNs, \(i) {i$target})))
# 403
allTF <- intersect(TF.df.co$Factor,
  union(allTF.source, allGene.target))

# * get sparseMatrix per subclass
GRNmats <- lapply(GRNs, function(grn) {
  i <- match(grn$source, allTF.source)
  j <- match(grn$target, allGene.target)
  x <- grn$coef
  sparseMatrix(i = i, j = j, x = x,
    dims = c(length(allTF.source), length(allGene.target)),
    dimnames = list(
      allTF.source, allGene.target), index1 = TRUE)
})
names(GRNmats) <- subclasses
saveRDS(GRNmats,
  file.path(outdir, "sa2.GRNmats.base_subclass.rds"))

# fit powerlaw
GRNmats <- readRDS(
  file.path(outdir, "sa2.GRNmats.base_subclass.rds"))
fitpList <- lapply(seq_along(GRNmats), function(i) {
  message("Powerlaw fitting: ", names(GRNmats)[i])
  fitPowerLaw(GRNmats[[i]])
})
names(fitpList) <- names(GRNmats)
bp.fitpower.df <- data.frame(
  subclass = names(fitpList),
  KSpvalue = vapply(fitpList, function(i){i$KS.p[1]}, FUN.VALUE = 1.0)
)

# * supple figure for powerlaw
points.text <- bp.fitpower.df[bp.fitpower.df$KSpvalue <= 0.05, ]
(
  bp.fitpwer <- ggplot(data = bp.fitpower.df, aes(x = 1, y = KSpvalue))  +
    geom_boxplot(outlier.colour="black", outlier.shape=16,
      outlier.size=2, notch=FALSE) +
    ggrepel::geom_text_repel(data = points.text,
      mapping = aes(x = 1 + 0.02, y = KSpvalue, label = subclass,
        color = "red"),
      size = 3) +
    theme_minimal() + 
    theme(
      axis.text.x = element_blank(),
      axis.title = element_text(colour = "black", size = 13),
      axis.text.y = element_text(colour = "black", size = 12),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = "NA"),
      panel.spacing = unit(0.002, "lines"),
      legend.position = "none"
    ) +
    ylab("p-value from Kolmogorov-Smirnov test")
)
ggsave(bp.fitpwer,
  filename = file.path(outdir, "Power-lab.distribution.fitting.pdf"), width = 3,
  height = 6)

# * eval how many TFs per gene, and how many genes per TF
nTFpGene <- lapply(GRNmats, function(i) {
  s <- as(i, "TsparseMatrix")
  s@x <- rep(1, length(s@i))
  r <- colSums(s)
  r[r!= 0]
})

nTFpGene <- unlist(nTFpGene)
boxplot(nTFpGene, outline = FALSE)

nGenepTF <- lapply(GRNmats, function(i) {
  s <- as(i, "TsparseMatrix")
  s@x <- rep(1, length(s@i))
  r <- rowSums(s)
  r[r!= 0]
})
nGenepTF <- unlist(nGenepTF)
boxplot(nGenepTF, outline = FALSE)

# * eval the size the networks
nTFpSubclass <- vapply(GRNmats, function(i) {
  s <- as(i, "TsparseMatrix")
  length(unique(s@i))
}, FUN.VALUE = 10)
boxplot(nTFpSubclass, outline = FALSE)

nGenepSubclass <- vapply(GRNmats, function(i) {
  s <- as(i, "TsparseMatrix")
  length(unique(s@j))
}, FUN.VALUE = 10)
boxplot(nGenepSubclass, outline = FALSE)



# * get TF2TF (TF2 for short)
TFs <- allTF.source
Genes <- allGene.target
## GRNBigTable <- t(vapply(GRNmats, function(grn) {
##   as.vector(grn)
## }, FUN.VALUE = rep(0.0, length(TFs) * length(Genes))))

## aGRNmat <- GRNmats[[1]]
## TFGenePairs <- vapply(colnames(aGRNmat), function(i) {
##   paste(rownames(aGRNmat), i, sep = "@")
## }, FUN.VALUE = rep("TF@Gene", length(TFs)))
## TFGenePairs <- as.vector(TFGenePairs)
## rownames(GRNBigTable) <- names(GRNs)
## colnames(GRNBigTable) <- TFGenePairs
## zeroIndex <- colSums(GRNBigTable) == 0
## GRNBigTable2 <- GRNBigTable[, !zeroIndex]

# * create TF2s
TF2s <- lapply(GRNs, toTF2, TFs = allTF.source)
names(TF2s) <- names(GRNs)
saveRDS(TF2s, file.path(outdir, "sa2.TF2s.list.rds"))

TF2s <- readRDS(file.path(outdir, "sa2.TF2s.list.rds"))
# * analyze TF2 network motif
TF2dir <- outdir
subclasses <- names(TF2s)
motifFile <- file.path(projdir, "05.GRN",
  "supple.07.02.list.network.motif.R")
# load defination of motfis
source(motifFile)
tag <- "color3"

library(future.apply)
plan(multicore, workers = 4)
options(future.globals.maxSize = 8e9)
cntMotifs.all <- future.apply::future_vapply(subclasses, function(s) {
  message("Working on subclass: ", s)
  TF2 <- TF2s[[s]]
  ## TF2 <- read.table(file = file.path(args$TF2dir,
  ##   paste0(s, ".TF2.csv")), header = FALSE, sep = ",",
  ##   quote = "")
  ## colnames(TF2) <- c("source", "target", "coef")
  g <- toiGraph3(TF2)
  cntMotifs <- countMotifs(netMotifList = netMotifList, g = g)
  motifs <- getMotifs(netMotifList = netMotifList, g = g)
  saveRDS(object = motifs, file = file.path(outdir,
    paste("sa2", s, "motifs", tag, "rds", sep = ".")))
  return(cntMotifs)
}, FUN.VALUE = rep(0, length(netMotifList)))

rownames(cntMotifs.all) <- names(netMotifList)
colnames(cntMotifs.all) <- subclasses
saveRDS(cntMotifs.all, file = file.path(outdir,
  paste0("sa2.all.subclass.cnt.motif.", tag, ".rds")))

cntMotifs.all <- readRDS(file.path(outdir,
  paste0("sa2.all.subclass.cnt.motif.", tag, ".rds")))

# * set subclass - main regions
scnmMap <- cembav2env$load.sa2.subclass.names()
rownames(scnmMap) <- scnmMap$subclass_id_label
sc2r <- read.table(
  file = file.path(projdir, "meta", "sa2.subclass2region2score.csv"),
  header = TRUE, sep = ",")
r2color <- S4Vectors::second(largeRegionColors)
names(r2color) <- S4Vectors::first(largeRegionColors)
sc2r$rcolor <- r2color[sc2r$r]
rownames(sc2r)

# * set subclass - nt type
atacMeta <- readRDS(cembav2env$sa2metaFile)
sc2nt <- unique(atacMeta[, c("subclass_label_v3", "NT_v3", "nt_type_color_v3")])
sc2nt$nt_type_color_v3[is.na(sc2nt$nt_type_color_v3)] <- "grey"
rownames(scnmMap) <- scnmMap$subclass_label
sc2nt$subclass <- scnmMap[sc2nt$subclass_label_v3, "subclass_label_peak"]
rownames(sc2nt) <- sc2nt$subclass
ntcolor <- sc2nt[, c("NT_v3", "nt_type_color_v3")] |>
  unique() |>
  x => {
    r <- x[, 2] |>
      setNames(object = _, nm = x[ ,1])
  }


# * set subclass order
sc2id <- scnmMap[, c("subclass_label_peak", "subclass_id")]
rownames(sc2id) <- sc2id$subclass_label_peak
atacsc_id <- sc2id[rownames(sc2nt), "subclass_id"]
osc <- rownames(sc2nt)[order(atacsc_id)]

# * prepare data format for ggplot
cnt <- readRDS(file.path(outdir,
  "sa2.all.subclass.cnt.motif.color3.rds"))

motri <- df2tri(as.data.frame(t(cnt[, osc[osc %in% colnames(cnt)]])))
colnames(motri) <- c("subclass", "motif", "count")
motri$mainclass <- sc2nt[motri$subclass, "NT_v3"]
motri$region <- sc2r[motri$subclass, "r"]
motri$regionpct <- sc2r[motri$subclass, "r_s"]

# histgram for main class
# can filter Chol (1), Sero(1), Glut-GABA (1), Dopa (2)
motri.class <- motri |>
  subset(x = _, subset = !(mainclass %in% c("Chol", "Sero", "Dopa", "Glut-GABA"))) |>
  subset(x = _, subset = !motif %in% c("NFL"))
(
  hist.mc <- ggplot(data = motri.class, aes(x = count, fill = mainclass)) +
  stat_density(aes(x = count, color = mainclass), na.rm = TRUE,
    show.legend = FALSE, geom = "line", position = "identity") +
  scale_color_manual(values = ntcolor) +
  geom_histogram(alpha = 0.6, position = "identity",
    na.rm = TRUE, bins = 30,
    aes(y = ..density..)) +
  scale_fill_manual(values = ntcolor) +
  facet_wrap(~motif, scales="free") +
   theme_bw()
)
ggsave(hist.mc, filename = file.path(outfigdir, "sa2.hist.motif.class.pdf"),
  width = 14, height = 10)
  
## case by case plotting
motifs <- unique(motri$motif)
outfigdir <- file.path(projdir, "figures", "out/fig4")
if(!dir.exists(outfigdir)) {
  dir.create(outfigdir)
}
for(i in seq_along(motifs)) {
  motri.tmp <- subset(motri, motif %in% motifs[i])
  hist.tmp <- ggplot(data =motri.tmp, aes(x = count, fill = mainclass)) +
    stat_density(aes(x = count, color = mainclass), na.rm = TRUE,
      show.legend = FALSE, geom = "line", position = "identity",
      linewidth = 0.5) +
    scale_color_manual(values = ntcolor)+
    geom_histogram(alpha = 0.4, position = "identity", na.rm = TRUE,
      bins = 25, aes(y = ..density..)) +
    scale_fill_manual(values = ntcolor) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
    ggtitle(motifs[i]) + 
    ## coord_cartesian() + 
    theme_bw() +
    theme(
      axis.title = element_text(colour = "black", size = 12),
      axis.text = element_text(colour = "black", size = 9),
      legend.position = "bottom",
      plot.title = element_text(colour = "black", size = 12)
    )
  hist.tmp
  ggsave(hist.tmp,
    filename = file.path(outfigdir,
      paste0(motifs[i], ".mainclass.pdf")
    ), width = 5, height= 5)
}

## histogram for regions
motr.region <- motri
motr.region <- motr.region[motr.region$regionpct >= 0.6, ]
motr.region <- motr.region[!motr.region$motif %in% "NFL", ]

(
  hist.r <- ggplot(data = motr.region, aes(x = count, fill = region)) +
    stat_density(aes(x = count, color = region), na.rm = TRUE,
      show.legend = FALSE, geom = "line", position = "identity") +
    scale_color_manual(values = r2color) +
    geom_histogram(alpha = 0.6, position = "identity", na.rm = TRUE,
      bins = 30, aes(y = ..density..)) +
    scale_fill_manual(values = r2color) +
    facet_wrap(~motif, scales="free") +
    theme_bw()
)
ggsave(hist.r, filename = file.path(outfigdir, "sa2.hist.motif.region.pdf"),
  width = 14, height = 10)

# case by case plotting
motifs <- unique(motri$motif)
for (i in 1:length(motifs)) {
motri.tmp <- subset(motri, motif %in% motifs[i])
hist.tmp <- ggplot(data =motri.tmp, aes(x = count, fill = region)) +
  stat_density(aes(x = count, color = region), na.rm = TRUE,
    show.legend = FALSE, geom = "line", position = "identity",
    linewidth = 0.5) +
  scale_color_manual(values = r2color)+
  geom_histogram(alpha = 0.4, position = "identity", na.rm = TRUE,
    bins = 25, aes(y = ..density..)) +
  scale_fill_manual(values = r2color) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  ggtitle(motifs[i]) + 
  ## coord_cartesian() + 
  theme_bw() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 9),
    legend.position = "right",
    plot.title = element_text(colour = "black", size = 12)
  )
hist.tmp
ggsave(hist.tmp,
  filename = file.path(outfigdir,
    paste0(motifs[i], ".region.pdf")
    ), width = 5, height= 4)
}

# ** heatmap on subclass with meta info (major region, size, class)
motif2color <- ArchRPalettes$calm[1:nrow(cnt)]
names(motif2color) <- rownames(cnt)
removeMotif <- c("SIM", "NFL")
sc2cnt <- t(cnt) |>
  as.data.frame() |>
  x => x[osc[osc %in% colnames(cnt)], ]
m2sc <- sc2cnt[, !colnames(sc2cnt) %in% removeMotif]
m2sc <- t(as.matrix(m2sc))
scnms <- colnames(m2sc)
m2sc <- m2sc %*% diag(1 / colSums(m2sc))
colnames(m2sc) <- scnms

motri.sc <- df2tri(df = m2sc)
colnames(motri.sc) <- c("motif", "subclass", "count")
# tmp remove MSN D1 since it has an empty bar at top
## motri.sc <- subset(motri.sc, !subclass %in% "MSN_D1_Gaba")
colnames(motri.sc) <- c("motif", "subclass", "count")
motri.sc$mc <- sc2nt[motri.sc$subclass, "NT_v3"]
motri.sc$subclass <- factor(motri.sc$subclass,
  levels = osc[osc %in% colnames(cnt)])

(
  bp.sc <- ggplot(data = motri.sc, aes(x = subclass, y = count, fill = motif)) +
  geom_bar(position = "stack", stat = "identity", na.rm = TRUE) +
  scale_fill_manual("Motif", values = motif2color) +
  scale_y_continuous(expand = expansion(0)) +
  ## scale_x_continuous(expansion(0)) +
  labs(x = "Cell subclass", y = "Ratio of motifs per cell subclass") +
  guides(fill = guide_legend(nrow = 2)) +
  ## ggnewscale::new_scale_fill() +
  ## geom_text(aes(label = "", x = subclass, y =-0.02), ) +
  ## geom_label(aes(fill = class, label = "", x = subclass, y =-0.02), size = 0.5) +
  ## scale_fill_manual("Class", values = class2color) + 
  ## scale_fill_manual("Motif and Class", values = c(motif2color, class2color)) +
  ## scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_classic()+
  theme(
    legend.position = "bottom",
    ## axis.title.x = element_blank(),
    ## axis.text.x = element_text(colour = "black", size = 6),
    axis.text.x = element_text(colour = "black",
      size = 6, angle = 90, hjust =0.95, vjust= 0.2),
    ## axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.title = element_text(colour = "black", size = 17),
    axis.ticks = element_blank()
  )
)
ggsave(bp.sc, filename = file.path(outfigdir,
  "motif.barplot.all.subclass.pdf"), width = 20, height = 8)

hist.motif <- ggplot(data = motri.sc, aes(x = count, fill = motif)) +
  geom_histogram(alpha = 0.5, position = "identity",
    na.rm = TRUE, bins = 60) +
  scale_fill_manual(values = motif2color) +
  theme_classic()
hist.motif

# * select motifs as examples
m <- "Fully_connected_triad"
motri[motri$motif == m, ] |> arrange(desc(count)) |> head(n=10)

motifFindingDir <- outdir
TF2Dir <- "TF2s.v1.nlogp"
sc <- "Peri_NN"
motifs <- readRDS(
  file.path(motifFindingDir, paste0("sa2.", sc,".motifs.color3.rds")))
unique(t(apply(motifs$Fully_connected_triad, 1, sort)))


sc <- "Astro-NT_NN"
motifs <- readRDS(
  file.path(motifFindingDir, paste0(sc,".motif.color3.rds")))
unique(t(apply(motifs$Fully_connected_triad, 1, sort)))

sc <- "VLMC_NN"
motifs <- readRDS(
  file.path(motifFindingDir, paste0(sc,".motif.color3.rds")))
unique(t(apply(motifs$Fully_connected_triad, 1, sort)))

sc <- "BAM_NN"
motifs <- readRDS(
  file.path(motifFindingDir, paste0("sa2.", sc,".motifs.color3.rds")))
unique(t(apply(motifs$Regulated_double_positive, 1, sort)))

sc <- "Endo_NN"
motifs <- readRDS(
  file.path(motifFindingDir, paste0(sc,".motif.color3.rds")))
unique(t(apply(motifs$Regulated_double_positive, 1, sort)))

sc <- "Ependymal_NN"
motifs <- readRDS(
  file.path(motifFindingDir, paste0(sc,".motif.color3.rds")))
unique(t(apply(motifs$Regulated_double_positive, 1, sort)))



# * check subclass having high numebr of motifs.
motifnm <- "Fully_connected_triad"
sort(cnt[motifnm, ], decreasing = TRUE)[1:5]


motifnm <- "Regulated_double_positive"
sort(cnt[motifnm, ], decreasing = TRUE)[1:5]


# * draw particular motif histogram
motifnm <- "Regulated_double_positive"
motri.tmp <-  motri.class[motri.class$motif %in% motifnm, ]
(
  hist.case.mc <- ggplot(data = motri.tmp, aes(x = count, fill = mainclass)) +
  stat_density(aes(x = count, color = mainclass), na.rm = TRUE,
    show.legend = FALSE, geom = "line", position = "identity") +
  scale_color_manual(values = ntcolor) +
  geom_histogram(alpha = 0.6, position = "identity",
    na.rm = TRUE, bins = 30,
    aes(y = ..density..)) +
  scale_fill_manual(values = ntcolor) +
   theme_bw()
)
ggsave(hist.case.mc, filename = file.path(outfigdir, "sa2.hist.regulated_double_positve.class.pdf"),
  width = 6, height = 4)

sc <- "BAM_NN"
motifs <- readRDS(
  file.path(motifFindingDir, paste0("sa2.", sc,".motifs.color3.rds")))
unique(t(apply(motifs$Regulated_double_positive, 1, sort)))

