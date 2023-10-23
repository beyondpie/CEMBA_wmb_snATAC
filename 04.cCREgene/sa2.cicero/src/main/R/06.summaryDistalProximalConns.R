library(ggplot2)
library(data.table)
library(cowplot)
library(scales)
packdir <- file.path(here::here(), "package", "R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env)
import::from(.from = "annot.R", .directory = packdir,
  getSubclassMeta.ordered, mapAllenSubclassTomCGSubclass2)
import::from(.from = "colors.R", .directory = packdir,
  ArchRPalettes)

# * themes
theme_my_minimal <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    axis.title = element_text(colour = "black", family = "serif", size = 10),
    axis.text = element_text(colour = "black", family = "serif", size = 10),
    plot.title = element_text(colour = "black", family = "serif",
      hjust = 0.5)
  )


# * configs
## load L3 meta
IntL3Meta <- data.table::fread(file = cembav2env$IntL3MetaFile,
  header = TRUE, data.table = FALSE)
rownames(IntL3Meta) <- IntL3Meta$L3
## order of subclass
L3.Intv2.hc <- readRDS(
  file.path(here::here(), "figures", "L3.hc.meta.rds"))
L3.hc <- L3.Intv2.hc[[4]]
order.L3 <- with(L3.hc, hc$labels[hc$order])
order.Allen.L3 <- IntL3Meta[order.L3, "L3AllenAnnot"]
order.Allen <- unique(order.Allen.L3)
order.AllenConcat <- mapAllenSubclassTomCGSubclass2(order.Allen)

## summary of conns
metaCol <- "AllenAnnotConcat"
fromdir <- file.path("out", metaCol)
allConnsStatFile <- file.path(fromdir,
  paste0("mba.whole.", metaCol, ".merge.conns.sta.all"))

allProxDistFile <- file.path(fromdir,
  paste0("mba.whole.",metaCol, ".merge.pdc.all"))
allProxDistStatFile <- file.path(fromdir,
  paste0("mba.whole.", metaCol, ".merge.pdc.sta.all"))

allGeneStatFile <- file.path(fromdir,
  paste0("mba.whole.", metaCol, ".merge.pdc.peak2gene.all"))
allPeakStatFile <- file.path(fromdir,
  paste0("mba.whole.", metaCol, ".merge.pdc.gene2peak.all"))
allDistFile <- file.path(fromdir,
  paste0("mba.whole.", metaCol, ".merge.pdc.dist.all"))
allPairFile <- file.path(fromdir,
  paste0("mba.whole.", metaCol, ".merge.bedpe.all"))

# * load data
conns.sta <- data.table::fread(
  file = allConnsStatFile, header = TRUE, data.table = FALSE, sep = "\t")
pdc.sta <- data.table::fread(
  file = allProxDistStatFile, header = TRUE, data.table = FALSE)
pdc.all <- data.table::fread(
  file = allProxDistFile, header = TRUE, data.table = FALSE)
pdc.dist.all <- data.table::fread(
  file = allDistFile, header = TRUE, data.table = FALSE)


# * distal-proximal, distal-distal, proximal-proximal ratio
boxplot.overall.df <- conns.sta[, 2:ncol(conns.sta)]
colors <- ArchRPalettes[["stallion"]][
  1:ncol(boxplot.overall.df)]
names(colors) <- colnames(boxplot.overall.df)
plist <- lapply(colnames(boxplot.overall.df), function(i) {
  t <- data.frame(
    name = rep(i, nrow(boxplot.overall.df)),
    count = boxplot.overall.df[, i])
  ylim1 <- boxplot.stats(t$count)$stats[c(1, 5)]
  ggplot(t, aes(x = name, y = count)) +
    geom_boxplot(fill = colors[i]) +
    labs(x = i) +
    coord_cartesian(ylim = ylim1*1.05) + 
    scale_y_continuous(label = scales::comma) +
    theme(
      axis.title = element_text(colour = "black", family = "serif", size = 15),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour = "black", family = "serif", size = 12),
      plot.title = element_text(colour = "black", family = "serif",
        hjust = 0.5)
    )
})
(p.bp.all <- cowplot::plot_grid(plotlist = plist, nrow = 1, labels = "auto"))

# * boxplot to show # peaks / gene

# * boxplot to show each subclass: #conns, #peaks, #genes in pdc
boxplot.overall.df <- pdc.sta[, 2:ncol(pdc.sta)]
colors <- ArchRPalettes[["stallion"]][
  1:ncol(boxplot.overall.df)]
names(colors) <- colnames(boxplot.overall.df)
plist <- lapply(colnames(boxplot.overall.df), function(i) {
  t <- data.frame(
    name = rep(i, nrow(boxplot.overall.df)),
    count = boxplot.overall.df[, i])
  ylim1 <- boxplot.stats(t$count)$stats[c(1, 5)]
  ggplot(t, aes(x = name, y = count)) +
    geom_boxplot(fill = colors[i]) +
    labs(x = i) +
    coord_cartesian(ylim = ylim1*1.05) + 
    scale_y_continuous(label = scales::comma) +
    theme(
      axis.title = element_text(colour = "black", family = "serif", size = 15),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour = "black", family = "serif", size = 12),
      plot.title = element_text(colour = "black", family = "serif",
        hjust = 0.5)
    )
})
(p.bp.all <- cowplot::plot_grid(plotlist = plist, nrow = 1, labels = "auto", ))
title <- cowplot::ggdraw() + draw_label("Proximal-distal conns", fontface = 'bold')
(p <- cowplot::plot_grid(title, p.bp.all, ncol = 1, rel_heights = c(0.1, 1)))

# * hist to show distances
median.dist <- median(pdc.dist.all$distance)
p.dist <- ggplot(data = pdc.dist.all, aes(x = distance)) +
  geom_histogram(fill = "lightblue", bins = 100) +
  labs(x = "Distance between proximal-distal connections",
    y = "Number of cCREs") +
  geom_vline(aes(xintercept = median(distance)),col='darkblue', linewidth=1) +
  ## geom_text(aes(label=median.dist,y=0,x=median.dist),
  ##           vjust=-1,col='black',size=5) + 
  scale_y_continuous(label = scales::comma) +
  scale_x_continuous(label = scales::comma) +
    theme(
      axis.title = element_text(colour = "black", family = "serif", size = 15),
      axis.text.x = element_text(colour = "black", family = "serif", size = 12),
      axis.text.y = element_text(colour = "black", family = "serif", size = 12),
      plot.title = element_text(colour = "black", family = "serif",
        hjust = 0.5))
p.dist
# * check with cembav1
