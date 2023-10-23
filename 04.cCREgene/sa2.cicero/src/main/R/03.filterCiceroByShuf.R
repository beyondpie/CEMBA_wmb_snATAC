library(optparse)
library(Matrix)
library(fitdistrplus)
library(ggplot2)

projdir <- here::here()
sa2_subclasspeak_dir <- file.path(projdir,
  "supple.07.peakcalling.allinone",
  "sa2.subclassv3.final.peak.srt")
pack_dir <- file.path(projdir, "package/R")
import::from(.from = "utils.R", .directory = pack_dir,
  checkArgsExistOrStop, prepareOutdir)
import::from(.from = "cembav2env.R", .directory = pack_dir,
  Sa2PeakCalling)


op <- list(
  make_option(c("--metaCol"), type = "character",
    default = "AllenAnnotConcat"),
  make_option(c("--group"), type = "character"),
  make_option(c("--ciceroDir"), type = "character"),
  make_option(c("--outdir"), type = "character"),
  make_option(c("--debug"), type = "integer", default = 0)
)

args <- parse_args(OptionParser(option_list = op))
debug <- if (args$debug > 0) {
  message("Under debug mode.")
  TRUE
} else {
  message("Under normal mode.")
  FALSE
}

# * test config
if (debug) {
  args$group <- "TU-ARH_Otp_Six6_Gaba"
  args$outdir <- "test_cicero"
  args$ciceroDir <- file.path(
    here::here(), "04.cCREgene", "sa2.cicero", "out/smamde_result")
}

# * check args
checkArgsExistOrStop(args)
if (!dir.exists(args$ciceroDir)) {
  stop(args$ciceroDir, " does not exist.")
}

# * prepare input files
conns.file <- file.path(args$ciceroDir, "cicero",
  paste0(args$group, "_cicero_connections.Rds"))
if (!file.exists(conns.file)) {
  stop(conns.file, " does not exist.")
}
shuf.conns.file <- file.path(args$ciceroDir, "shuf_cicero",
  paste0(args$group, "_cicero_shuffle_connections.Rds"))

filterByShuf <- TRUE
if (!file.exists(shuf.conns.file)) {
  message(shuf.conns.file, " does not exist.")
  warning("No filtering on", args$group)
  filterByShuf <- FALSE
}

# * prepare outfiles
prepareOutdir(args$outdir)
outprefix <- file.path(args$outdir, paste(
  "cicero", args$metaCol, args$group, sep = "."
))
outFitParamFile <- paste(outprefix, "fitConns.para.txt", sep = ".")
outFitResFile <- paste(outprefix, "fitConns.res.txt", sep = ".")
outFitDistFile <- paste(outprefix, "fitConns.dist.pdf", sep = ".")

## if ((file.exists(outFitParamFile)
##   && file.exists(outFitResFile))
##   && file.exists(outFitDistFile)) {
##   warning("quite the filter job since files exist.")
##   message(outFitParamFile, outFitResFile,
##     outFitDistFile, " exist, and quit the job.")
##   quit(save = "no", status = 0)
## }

# * get group-related peak coords
bpmatFile <- file.path(sa2_subclasspeak_dir,
  paste0(args$group, ".bed"))
if (!file.exists(bpmatFile)) {
  stop(bpmatFile, "does not exist.")
}
peaks <- Sa2PeakCalling$readBed4File(fnm = bpmatFile)
peakCoordsInGroup <- with(peaks,
  paste(chr, paste(start, end, sep = "_"), sep = "_"))

# * funcitons
filter1 <- function(raw.conns) {
  idx <- intersect(
    which(raw.conns$Peak1 %in% peakCoordsInGroup),
    which(raw.conns$Peak2 %in% peakCoordsInGroup))
  conns <- raw.conns[idx, ]
  # * dedup same pairs
  peak.pairs <- with(conns, paste(Peak1, Peak2, sep = "."))
  conns <- conns[!duplicated(peak.pairs), ]
  return(conns)
}

# * load co-accessible peaks
## raw.conns <- data.table::fread(
##   conns.file, sep = "\t", header = TRUE, data.table = FALSE)
raw.conns <- readRDS(conns.file)
raw.conns <- raw.conns[!is.na(raw.conns$coaccess), ]

# * filter based on peaks in the groups
conns <- filter1(raw.conns)

if (filterByShuf) {
  message("Filter cicero based on shuffled file for ", args$group)
  shuf.conns <- readRDS(shuf.conns.file)
  shuf.conns <- shuf.conns[!is.na(shuf.conns$coaccess), ]

  # * filter based on peaks in the groups
  shuf.conns <- filter1(shuf.conns)

  # * filter based on shuf signals
  fitnorm.shuf <- fitdistrplus::fitdist(
    data = shuf.conns$coaccess, distr = "norm", keepdata = FALSE)
  mean.shuf <- fitnorm.shuf$estimate[[1]]
  std.shuf <- fitnorm.shuf$estimate[[2]]

  # ** output fitting plots
  withr::with_pdf(new = outFitDistFile, code = {
    plot(fitnorm.shuf)
    hist(shuf.conns$coaccess, pch = 20, breaks = 25, prob = TRUE,
      main = paste(args$group, " shuffled co-accessibility scores"))
    dat <- rbind(data.frame(
      coaccess = shuf.conns$coaccess,
      class = rep("shuf", nrow(shuf.conns))
    ), data.frame(
      coaccess = conns$coaccess,
      class = rep("conns", nrow(conns))
    ))
    print(ggplot(dat, aes(x = coaccess, color = class)) +
            geom_density() +
            theme_bw() +
            #  xlim(-0.2, 0.2) +
            #  ylim(0,200) +
            theme(axis.text.x = element_text(colour = "black", size = 8),
              axis.text.y = element_text(colour = "black", size = 8),
              axis.title.x = element_text(colour = "black", size = 8),
              axis.title.y = element_text(colour = "black", size = 8)))
  })

  # ** save params
  para <- data.frame(
    group = args$group,
    metaCol = args$metaCol,
    meanShuf = as.numeric(mean.shuf),
    stdShuf = as.numeric(std.shuf)
  )
  write.table(para,
    outFitParamFile, quote = FALSE,
    col.names = TRUE, row.names = FALSE, sep = "\t")

  # ** get one-side pvalue and FDR for conns
  cal_p <- function(x, mean, sd, lower = TRUE) {
    p <- if (lower) {
      1 - pnorm(x, mean = mean, sd = sd, lower.tail = TRUE, log.p = FALSE)
    } else {
      1 - pnorm(x, mean = mean, sd = sd, log.p = FALSE)
    }
    return(p)
  }

  testPval <- vapply(conns$coaccess,
    cal_p, mean = mean.shuf, sd = std.shuf, lower = TRUE,
    FUN.VALUE = 0.0)
  testFDR <- p.adjust(testPval, method = "BH")
  conns$nlog10p <- (-log10(testPval))
  conns$FDR <- testFDR
  conns$p <- testPval
} else {
  conns$nlog10p <- 100000
  conns$FDR <- 0
  conns$p <- 0
}

write.table(
  conns, outFitResFile,
  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

message("Done")
