library(optparse)

packdir <- file.path(here::here(), "package/R")
import::from(
  .from = "utils.R",
  .directory = packdir,
  checkArgsExistOrStop, prepareOutdir,
  checkFileExistOrStop, colFn)
import::from(.from = "cembav2env.R",
  .directory = packdir, cembav2env)
import::from(
  .from = "hc.R",
  .directory = packdir,
  getNMFHeatMap)
import::from(
  .from = "annot.R",
  .directory = packdir,
  getSubclassMeta.ordered)
import::from(
  .from = "peak.R",
  .directory = packdir,
  loadStatPeak.NMF,
  processCPM.log2.scale.cap,
  downsample.index.1,
  downsampleStatPeak.1,
  orderModules.cpm.statPeak,
  orderStatPeakByModuleRank,
  getPeakModuleAvgScore.NMF)

op <- list(
  make_option(c("--nmfDir"), type = "character",
    default = "nmf_subclassIntv2_novlp_tscc/out"),
  make_option(c("--module"), type = "integer", default = 150),
  make_option(c("--tag"), type = "character", default = "novlp.Intv2")
)
args <- parse_args(OptionParser(option_list = op))
checkArgsExistOrStop(args)

if (!dir.exists(args$nmfDir)) {
  stop(args$nmfDir, " does not exist.")
}

# * configs
module <- args$module
tag <- args$tag

statPeakFile <- file.path(
  args$nmfDir,
  paste("nmfPmat", tag, paste0("r", module),
    "n0.statW", sep = "."))
if(!file.exists(statPeakFile)) {
  stop(statPeakFile, " does not exist.")
}

peaks <- data.table::fread(cembav2env$peakBedFile,
  header = FALSE, sep = "\t", data.table = FALSE)
colnames(peaks) <- c("chrom", "start", "end", "name")
peakCoords <- with(peaks,
  paste(chrom, paste(start, end, sep = "-"), sep = ":"))
names(peakCoords) <- peaks$name

cluster.order.hc <- read.table(cembav2env$subclass.order.hc.Intv2,
  header = FALSE, sep = ",")$V1
## TODO: update
ha_row <- NULL

# * load nmf result
statPeak <- loadStatPeak.NMF(file = statPeakFile)
n_total <- 10000
n_module <- length(unique(statPeak$moduleN))
n_avg <- floor(n_total / n_module)
statPeak.ds <- downsampleStatPeak.1(
  statPeak = statPeak,
  mod.ordered = NULL,
  size = n_avg,
  seed = 2022)

# * load cpm
cpm <- readRDS(cembav2env$subclassPmatCPMIntv2File)
cpm.log2 <- processCPM.log2.scale.cap(
  cpm = cpm,
  islog2 = TRUE,
  isScale = FALSE,
  isCapped = FALSE)
cpm.plot <- if (grepl("chr", statPeak.ds$peak[1])) {
  message("StatPeak coords start with chr.")
  cpm.log2[cluster.order.hc, statPeak.ds$peak]
}  else {
  message("StatPeak coords start with peak name.")
  cpm.log2[
    cluster.order.hc, peakCoords[statPeak.ds$peak]]
}

# * set module order
nmf.order <- orderModules.cpm.statPeak(
  cpm = cpm.plot,
  statPeak = statPeak.ds)

# * re-order statPeak
statPeak.ordered <- orderStatPeakByModuleRank(
  statPeak = statPeak.ds,
  mod.ordered = nmf.order)

# * reorder peaks by putting the global modules firstly
modules <- unique(statPeak.ordered$moduleN)
cpm.mod <- getPeakModuleAvgScore.NMF(
  cpm = cpm.plot,
  moduleOfPeak = statPeak.ordered$moduleN,
  fn = rowMeans
)
high.value <- quantile(cpm.plot, 0.95)
nsubclass.peak.plot <- colSums(cpm.mod >= high.value)
global_modules<- colnames(cpm.mod)[
  which(nsubclass.peak.plot >= (0.9 * nrow(cpm.plot)))]
mod.reorder <- c(
  global_modules, setdiff(colnames(cpm.mod), global_modules)
)
statPeak.reorder <- orderStatPeakByModuleRank(
  statPeak.ordered, mod.reorder
)

cpm.plot.reorder <- if (grepl("chr", statPeak.ds$peak[1])) {
  cpm.plot[ , statPeak.reorder$peak]
}  else {
  cpm.plot[ , peakCoords[statPeak.reorder$peak]]
}

# * set heatmap column annotation
moduleColor <- grDevices::colorRampPalette(
  colors = c("bisque", "burlywood4"))(length(mod.reorder))
names(moduleColor) <- nmf.order
ha_col <- ComplexHeatmap::HeatmapAnnotation(
  Module = statPeak.reorder$moduleN,
  FeatureScore = ComplexHeatmap::anno_barplot(
    statPeak.reorder$featureScore),
  col = list(
    Module = moduleColor
  ),
  which = "column",
  show_legend = c(FALSE))

# * get heatmap
p.cpm.log2 <- getNMFHeatMap(
  cpm.plot = cpm.plot.reorder,
  ## NOTE: update later
  ha_row = NULL,
  ha_col = ha_col,
  fontsize = 6,
  low.val.col = quantile(cpm.plot, 0.01),
  high.val.col = quantile(cpm.plot, 0.99))

# * plot figures
withr::with_pdf(
  new = file.path(args$nmfDir,
    paste(tag, paste0("r", module), "pdf", sep = ".")),
  code = {
    print(p.cpm.log2)
  }, height = 28, width = 20
)
