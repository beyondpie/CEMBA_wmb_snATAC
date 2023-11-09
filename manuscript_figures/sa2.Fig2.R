source("paper.R")
# * configs
atacMeta <- readRDS(cembav2env$sa2metaFile)
peakBed <- Sa2PeakCalling$loadFinalPeakBed()
nds <- 10000

# * load atac and mCG signals.
logCPM.scbyp <- Sa2PeakCalling$read.sa2CPM.pbysc(
  cap = 0.99)
mCG.all <- readRDS(Sa2PeakCalling$sa2.mCG.rds.File)

scnmMap <- cembav2env$load.sa2.subclass.names()
rownames(scnmMap) <- scnmMap$subclass_label_mCG
new_scnms <- scnmMap[rownames(mCG.all), "subclass_label_peak"]
rownames(mCG.all) <- new_scnms
# 244 subclasses shared between mCG and ours

scs <- intersect(rownames(logCPM.scbyp), rownames(mCG.all))
# sort scs based on allen order
scs_id <- scnmMap$subclass_id[match(scs, scnmMap$subclass_label_peak)]
scs_ord <- scs[order(scs_id)]

mCG <- mCG.all[scs_ord, ]
# remove NA part
colNA <- vapply(seq_len(ncol(mCG)), function(i) {
  r <- if (any(is.na(mCG[, i]))) { TRUE } else { FALSE }
}, FUN.VALUE = FALSE)
mCG <- mCG[, !colNA]
# top 2 with lowest mCG signals.
# 095 DMH Prdm13 Gaba	| #cell: 845 | # cell 32
# 153 MG-POL-SGN Nts Glut	| #cell: 414 | #cell 34

## 32.00   175.75   431.50   830.25 21362.00 
## > sort(sc2n.mCG)[1:10]
## SubClass
##          DMH Prdm13 Gaba      MG-POL-SGN Nts Glut     CEA-BST Gal Avp Gaba 
##                       32                       34                       35 
##        MM-ant Foxb1 Glut   TMv-PMv Tbx3 Hist-Gaba MDRN Hoxb5 Ebf2 Gly-Gaba 
##                       35                       37                       38 
##             PAG Ucn Glut        NLL-SOC Spp1 Glut       DMX VII Tbx20 Chol 
##                       41                       44                       46 
## PRNc Prox1 Brs3 Gly-Gaba 
##                       46 

mCGMeta <- data.table::fread(
  file = file.path(projdir, "12.DMR", "src/main/resource",
    "CEMBA.mC.Metadata.csv"),
  sep = ",",
  header = TRUE,
  data.table = FALSE
)

sc2n.mCG <- with(mCGMeta, table(SubClass))


logCPM <- logCPM.scbyp[scs_ord, ]


# * Fig.2f top panel
# load nmf result
allStatPeak <- loadStatPeak.NMF(
  file.path(projdir, "figures", "data",
    "fig2/sa2_allpeak_nmf", "nmfPmat.allpeak.r150.n0.statW")
)
nmod <- length(unique(allStatPeak$moduleN))

## generate heatmap
allStatPeak.ds <- downsampleStatPeak.1(
  statPeak = allStatPeak,
  mod.ordered = NULL,
  size = floor(nds / nmod),
  seed = 2022
)
peak.ds <- intersect(allStatPeak.ds$peak, colnames(mCG))
mCG.ds <- mCG[, peak.ds]
logCPM.ds <- logCPM[, peak.ds]
allStatPeak.ds <- allStatPeak.ds[allStatPeak.ds$peak %in% peak.ds, ]

# get global modules
signalScByMod <- getPeakModuleAvgScore.NMF(
  cpm = logCPM.ds,
  moduleOfPeak = allStatPeak.ds$moduleN,
  fn = rowMeans)
global_modules <- getGlobalModule(
  signalScByMod,
  quantileHighSignal = 0.95,
  ratioOfSubclass = 0.9)

# set module order
mod.ord <- orderModules.cpm.statPeak2(
  cpm = logCPM.ds,
  statPeak = allStatPeak.ds)
allStatPeak.ord <- orderStatPeakByModuleRank(
  statPeak = allStatPeak.ds, mod.ordered = mod.ord)
modules <- unique(allStatPeak.ord$moduleN)

mod.ord <- c(global_modules,
  setdiff(modules, global_modules))
data.table::fwrite(
  x = data.frame(mod.order = mod.ord), file = file.path(
  projdir, "figures/data/fig2", "sa2.allpeak.nmf.module.order.txt"),
  col.names = FALSE)

statPeak.ord <- orderStatPeakByModuleRank(
  allStatPeak.ord, mod.ord)
# then reorderagain
logCPM.plot <- logCPM.ds[, statPeak.ord$peak]
mCG.plot <- mCG.ds[, statPeak.ord$peak]
## ComplexHeatmap::draw(heatmap.all.atac,
##   heatmap_legend_side = "bottom")

hm.all <- combineHeatmap(logCPM = logCPM.plot,
  mCG = mCG.plot)
withr::with_pdf(
  new = file.path(projdir, "figures", "out/figures",
    "nmf_all.pdf"),
  code = { print(hm.all) },
  height = 15, width = 20
)


# * Fig.2f bottom panel
nonOvlpBed <- Sa2PeakCalling$loadNonOvlpPeakBed()
novlpStatPeak <- loadStatPeak.NMF(
  file.path(projdir, "figures", "data",
    "fig2/sa2_novlpeak_nmf", "nmfPmat.novlpeak.r150.n0.statW")
)
nmod <- length(unique(novlpStatPeak$moduleN))

## generate heatmap
novlpStatPeak.ds <- downsampleStatPeak.1(
  statPeak = novlpStatPeak,
  mod.ordered = NULL,
  size = floor(nds / nmod),
  seed = 2022
)
peak.ds <- intersect(novlpStatPeak.ds$peak, colnames(mCG))
mCG.ds <- mCG[, peak.ds]
logCPM.ds <- logCPM[, peak.ds]
novlpStatPeak.ds <- novlpStatPeak.ds[novlpStatPeak.ds$peak %in% peak.ds, ]

# set module order
mod.ord <- orderModules.cpm.statPeak2(
  cpm = logCPM.ds,
  statPeak = novlpStatPeak.ds)
novlpStatPeak.ord <- orderStatPeakByModuleRank(
  statPeak = novlpStatPeak.ds, mod.ordered = mod.ord)
# then reorderagain
logCPM.plot <- logCPM.ds[, novlpStatPeak.ord$peak]
mCG.plot <- mCG.ds[, novlpStatPeak.ord$peak]
## ComplexHeatmap::draw(heatmap.all.atac,
##   heatmap_legend_side = "bottom")

hm.novlp <- combineHeatmap(logCPM = logCPM.plot,
  mCG = mCG.plot)
withr::with_pdf(
  new = file.path(projdir, "figures", "out/figures",
    "nmf_novlp.pdf"),
  code = { print(hm.novlp) },
  height = 15, width = 20
)

# * Fig2d
## construct pL4 binary pmat.
pL4s <- unique(atacMeta$pL4)
pL4peakdir <- file.path(projdir,
  "supple.07.peakcalling.allinone", "sa2.pL4.final.peak.srt")
pL4Peaks <- lapply(pL4s, \(pl4) {
  fnm <- file.path(pL4peakdir, paste0(pl4, ".bed"))
  Sa2PeakCalling$readBed4File(fnm)
})
names(pL4Peaks) <- pL4s

bpmat <- vapply(pL4s, \(x) {
  r <- rep(0, nrow(peakBed))
  r[match(pL4Peaks[[x]]$name , peakBed$name)] <- 1
  return(r)
}, rep(0, nrow(peakBed))) |> t()
rownames(bpmat) <- pL4s
colnames(bpmat) <- peakBed$name
saveRDS(bpmat,
  file.path(projdir, "figures", "data", "sa2.binary.pmat.pL4.rds")
)
bpmat <- readRDS(file.path(
  projdir, "figures", "data", "sa2.binary.pmat.pL4.rds"
))
noccurPeak <- colSums(bpmat >= 1)
noccurNonOvlp <- noccurPeak[nonOvlpBed$name]
noccurOvlp <- noccurPeak[!(names(noccurPeak) %in% nonOvlpBed$name)]

countpL4 <- data.frame(
  group = c(rep("Ovlp-DHS", length(noccurOvlp)),
    rep("NonOvlp-DHS", nrow(nonOvlpBed))),
  count = c(noccurOvlp, noccurNonOvlp)
)
getPieChartData <- function(count, maxVal = 9, class = "nonovlp") {
  s <- rle(sort(count))
  r <- s$lengths[s$values <= maxVal]
  a <- s$lengths[s$values > maxVal]
  r <- c(r, sum(a))
  names(r) <- c(as.character(1:maxVal), paste0(">=", maxVal + 1))
  df <- data.frame(count = names(r), num = r)
  df$percent <- round(df$num / sum(df$num), digits = 2) * 100
  df$ypos <- cumsum(df$percent) - 0.5 * df$percent
  df$count <- factor(df$count, levels = df$count)
  df$class <- rep(class, nrow(df))
  return(df)
}

maxVal <- 9
nonOvlpPC <- getPieChartData(
  count = countpL4$count[countpL4$group %in% "NonOvlp-DHS"],
  maxVal = maxVal,
  class = "nonovlp")
ovlpPC <- getPieChartData(
  count = countpL4$count[countpL4$group %in% "Ovlp-DHS"],
  maxVal = maxVal,
  class = "ovlp")
piedf <- rbind(ovlpPC, nonOvlpPC)
colors <- ArchRPalettes$stallion[as.character(1:10)]
names(colors) <- NULL
pieLabel <- c("nonovlp" = "NonOvlp", "ovlp" = "Ovlp")
(
  piechart <- ggplot(data = piedf,
    aes(x = "", y = percent, fill = count)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    guides(fill = guide_legend(title = "# of cell types")) +
    scale_fill_manual(values = colors) +
    facet_grid(cols = vars(class), labeller = labeller(class=pieLabel)) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 7, colour = "black"),
    legend.key.size = unit(4, 'mm'),
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
ggsave(filename = file.path(projdir, "figures",
  "out/figures", "sa2.pieChar_ncelltypes.pdf"),
  plot = piechart)

sum(noccurOvlp >= 10) / length(noccurOvlp)
sum(noccurOvlp == 1) / length(noccurOvlp)
sum(noccurNonOvlp >= 10) / length(noccurNonOvlp)
sum(noccurNonOvlp == 2) / length(noccurNonOvlp)
sum(noccurNonOvlp == 3) / length(noccurNonOvlp)


# * Others
# ** Update sa2.subclass.names with mCG names
## scnmMap$subclass_label_mCG <- gsub(
##   " ", "_", scnmMap$subclass_label) |>
##   gsub("/", "_", x = _)
## write.table(scnmMap, file = cembav2env$sa2.subclass.names,
##   col.names = TRUE, sep = ",", quote = FALSE, row.names = FALSE)

# ** output global_modules'related peaks
## peaksGlobal <- novlpStatPeak$peak[
##   novlpStatPeak$moduleN %in% global_modules]
## peaksdf <- vapply(peaksGlobal, str_split_1,
##   FUN.VALUE = rep("1",3), pattern = ":|-") |>
##   t() |>
##   as.data.frame()
## colnames(peaksdf) <- c("chr", "start", "end")
## peaksdf$name <- peaksGlobal
## write.table(peaksdf,
##   file = file.path(projdir, "20.nmf", "out",
##     "peaks.from.global.modules.bed"), sep = "\t",
##   col.names = FALSE, row.names = FALSE, quote = FALSE)


# * peaks per cell types
nL4s <- table(atacMeta$pL4)
bpmat <- readRDS(
  file.path(projdir, "figures", "data", "sa2.binary.pmat.pL4.rds")
)
npL4s <- rowSums(bpmat)

n2pL4 <- data.frame(
  nL4 = as.vector(nL4s),
  npL4 = npL4s[names(nL4s)]
)


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

(p.n2p <- ggplot(data = n2pL4, aes(x = nL4, y = npL4)) +
   geom_point(size = 2, shape = 19) +
   coord_trans(x="log10", y="log10") +
   scale_y_continuous(label=scientific_10, trans = 'log10', n.breaks = 10) + 
   scale_x_continuous(trans='log10', n.breaks = 10, label=scientific_10) +
   xlab("# of nuclei in L4-level cluster") +
   ylab("# of peak")+
   theme_bw() +
   theme(
     axis.text.y = element_text(colour = "black", size = 14),
     axis.text.x = element_text(colour = "black", size = 13, angle = 40, hjust=1),
     axis.title = element_text(colour = "black", size = 15)
   )
)
ggsave(p.n2p, filename = file.path(
  outfigdir, "sa2.npeak.vs.size.pL4.pdf"
), width = 10, height = 5)

# * plot accessibility's median vs variation
outfigdir <- file.path(projdir, "figures", "out",
  "figure2")
if(!dir.exists(outfigdir)) {
  dir.create(outfigdir)
}
nonOvlpPeaknm <- nonOvlpBed$name
ovlpPeaknm <- setdiff(colnames(logCPM.scbyp), nonOvlpPeaknm)
delta.pmat <- apply(logCPM.scbyp, 2, max) - apply(logCPM.scbyp, 2, min)
median.pmat <- apply(logCPM.scbyp, 2, median)

delta.pmat.dhs <- delta.pmat[ovlpPeaknm]
delta.pmat.ndhs <- delta.pmat[nonOvlpPeaknm]

median.pmat.dhs <- median.pmat[ovlpPeaknm]
median.pmat.ndhs <- median.pmat[nonOvlpPeaknm]

buylrd <- c("#313695", "#4575B4", "#74ADD1",
  "#ABD9E9", "#E0F3F8", "#FFFFBF",
  "#FEE090", "#FDAE61", "#F46D43",
  "#D73027", "#A50026")

myColRamp <- colorRampPalette(c(buylrd))

tag <- "subclass"
withr::with_pdf(
  new = file.path(outfigdir,
    "peak.minmax.scatter.ovlp2novlp.pdf"), code = {
      withr::with_par(new = list(
        mar = c(5, 4, 4, 5) + .1,
        mfrow = c(1, 2),
        cex = 2
      ),
      code = {
        smoothScatter(
          x = delta.pmat.dhs, y = median.pmat.dhs,
          colramp = myColRamp,
          main = paste("Accessibility of DHS-ovlp cCREs on\n",
            tag, "level", sep = " "),
          xlab = latex2exp::TeX(input = "$\\Delta$ accessibility"),
          ylab = "Median accessibility",
          ## postPlotHook = fudgeit,
          nrpoints = 0,
          nbin = 256)
        smoothScatter(
          x = delta.pmat.ndhs, y = median.pmat.ndhs,
          colramp = myColRamp,
          main = paste("Accessibility of non-DHS-ovlp cCREs on\n",
            tag, "level", sep = " "),
          xlab = latex2exp::TeX(input = "$\\Delta$ accessibility"),
          ylab = "Median accessibility",
          ## postPlotHook = fudgeit,
          nrpoints = 0,
          nbin = 256)
      })
    },
  width = 20, height = 10)
