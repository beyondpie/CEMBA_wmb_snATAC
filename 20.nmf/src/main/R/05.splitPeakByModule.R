library(data.table)
library(optparse)
packdir <- file.path(here::here(), "package/R")
import::from(.from = "utils.R", .directory = packdir,
  checkArgsExistOrStop, prepareOutdir, checkFileExistOrStop)
import::from(.from = "peak.R", .directory = packdir,
  loadStatPeak.NMF)
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env)

op <- list(
  make_option(c("--nmfDir"), type = "character",
    default = "nmf_ppdc/out"),
  make_option(c("--module"), type = "integer",
    default = 54),
  make_option(c("--tag"), type = "character", default = "ppdc")
)

args <- parse_args(OptionParser(option_list = op))
checkArgsExistOrStop(args)

if(!dir.exists(args$nmfDir)) {
  stop(args$nmfDir, " does not exist.")
}

mod.nmf <- args$module
tag <- args$tag

outDir <- file.path(args$nmfDir,
  paste("nmf", tag, paste0("r", mod.nmf), "motif", sep = "."))
prepareOutdir(outDir)

# * functions
convertPeakToBed <- function(peakBed, peaknms, outFile = NULL) {
  r <- peakBed[peaknms, ]
  if(!is.null(outFile)) {
    write.table(x = r, file = outFile, quote = FALSE, sep = "\t",
      row.names = FALSE, col.names = FALSE)
  }
  return(r)
}

# * load peaks
peakBed <- data.table::fread(cembav2env$peakBedFile,
  header = FALSE, sep = "\t", data.table = FALSE)
colnames(peakBed) <- c("chrom", "start", "end", "name")
rownames(peakBed) <- peakBed$name

# * nmf modules 
nmfPeakStat <- loadStatPeak.NMF(
  file = file.path(args$nmfDir,
    paste("nmfPmat", tag,
      paste0("r", mod.nmf), "n0", "statW", sep = ".")))

modules <- unique(nmfPeakStat$class0 + 1)

# * save peaks from each module to a seperate bed file
invisible(lapply(modules, function(i) {
  outFile <- file.path(outDir,
    paste0("r", mod.nmf, "_n", i, ".cCREs.bed"))
  message("Writing peak bed file to: ", outFile)
  peaks <- with(nmfPeakStat, peak[class0 == (i-1)])
  convertPeakToBed(peakBed = peakBed, peaknms = peaks, outFile = outFile)
}))





