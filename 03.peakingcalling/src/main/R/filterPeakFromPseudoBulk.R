library(tidyverse)
library(rlang)
library(R6)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

library(data.table)
library(GenomicRanges)
library(S4Vectors)
library(BSgenome)


rdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env)
import::from(.from = "peak.R", .directory = rdir,
  mapL4pc2L4)

work_dir <- file.path(here::here(), "18.snap2_peakcalling")
outdir4check <- file.path(work_dir, "out/check_pc")

# soft threshold
## we choose a cutoff of FDR less than 0.1, 0.05, 0.025, 0.01, and 0.001,
## corresponding to the situations when the number of reads is with then
## range of 0-5 million, 5-25 million, 25-50 million, 50-100 million, and
## 100 million and above.
# for each pL4, let's get the softthreshold vs SPM threshold

# * load union peaks
unionPeaks <- data.table::fread(
  file = file.path(outdir4check, "mba.whole.filteredNfixed.union.peakSet"),
  header = TRUE, sep = "\t",
  data.table = FALSE
)

# * load meta
atacMeta <- readRDS(cembav2env$sa2metaFile)

# * check all the L4 cls have been covered by pc.
uL4s <- unique(atacMeta$L4)
upL4s <- unique(unionPeaks$label)
pL4_L4 <- mapL4pc2L4(upL4s)
mappedL4 <- pL4_L4$L4 |> gsub(pattern = "^.*\\.", "", x = _)
pL4_L4$L4 <- mappedL4
rownames(pL4_L4) <- pL4_L4$L4
uL4s[!(uL4s %in% mappedL4)]

# * check how many l4 cls / pcls are not balanced
neuL4pcMeta <- read.csv(
  file = file.path(work_dir,
    "src/main/resource", "neuron_L4pc2sizes_cca-k50-cl_v1.csv"),
  header = FALSE, sep = ",")
nnL4pcMeta <- read.csv(
  file = file.path(work_dir,
    "src/main/resource", "nn_L4pc2sizes_cca-k49-cl_v1.csv"),
  header = FALSE, sep = ",")

L4pcMeta <- rbind(neuL4pcMeta, nnL4pcMeta) |>
  setNames(object = _, nm = c("pL4", "size", "n1", "n2"))
rownames(L4pcMeta) <- L4pcMeta$pL4
L4pcMeta$r1 <- with(L4pcMeta, n1 / size)
L4pcMeta$r2 <- with(L4pcMeta, n2 / size)
L4pcMeta[(L4pcMeta$r1 < 0.2) | (L4pcMeta$r2 < 0.2), ]

# * get size and nfrag per pL4
b2pL4 <- data.frame(
  barcode = atacMeta$barcode2,
  pL4 = pL4_L4[atacMeta$L4, "pL4"],
  nfrag = atacMeta$n_fragment
)
## t <- as.data.frame(table(b2pL4$pL4))
## pL4Meta <- data.frame(
##   pL4 = gsub("^.*\\.", "", as.character(t$Var1)),
##   size = t$Freq
## )
## rownames(pL4Meta) <- pL4Meta$pL4

pL4nread <- b2pL4 |> group_by(pL4) |>
  summarise(rmillion = sum(nfrag) / 1e6) |>
  as.data.frame()
# gL4 means grouped L4, like nn.[cl] or neuron.[cl]
# during peak calling
pL4nread$gL4 <- pL4nread$pL4
pL4nread$pL4 <- gsub("^.*\\.", "", pL4nread$pL4)
rownames(pL4nread) <- pL4nread$pL4
L4pcMeta$rmillion <- pL4nread[L4pcMeta$pL4, "rmillion"]
L4pcMeta$gL4 <- pL4nread[L4pcMeta$pL4, "gL4"]
write.table(L4pcMeta,
  file = file.path(work_dir, "src/main/resource",
    "all.pL4.meta.csv"),
  quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

# * check SPM and soft threshold
unionPeaks2 <- filterBySPM(unionPeaks, 2)
unionPeaks5 <- filterBySPM(unionPeaks, 5)

# * filter either SPM or softthreshold per cluster.
tostr.PeakSet <- function(peakset, sep1 = ":", sep2 = "-") {
  paste(peakset[, 1],
    paste(peakset[, 2], peakset[, 3], sep = sep2),
    sep = sep1)
}

toGR.PeakSet <- function(peakset) {
  if (is.null(colnames(peakset))) {
    stop("set peakSet colnames before use it.")
  }
  gr <- GenomicRanges::GRanges(
    seqnames = peakset[, 1],
    ranges = IRanges::IRanges(
      start = peakset[, 2],
      end = peakset[, 3]
  ))
  if (ncol(peakset) > 3) {
    for (j in seq(from = 4, to = ncol(peakset))) {
      S4Vectors::mcols(gr)[, colnames(peakset)[j]] <-
        peakset[, j]
    }
  }
  return(gr)
}

#' @param nfragm float nfragment in million
getLog10NegQvalThreshold <- function(nfragm) {
  if (nfragm < 5) {
    return(-log10(0.1))
  } else if (nfragm < 25) {
    return(-log10(0.05))
  } else if (nfragm < 50) {
    return(-log10(0.025))
  } else if (nfragm < 100) {
    return(-log10(0.01))
  } else {
    return(-log10(0.001))
  }
}

L4PeakSet <- R6::R6Class(classname = "L4PeakSet",
  public = list(
    peakset = data.frame(),
    ncell = 0L,
    nfragm = 0.1,
    minMacs2Score = -log10(0.01),
    minSPM = 2.0,
    peaksetSoftfilter = data.frame(),
    peaksetSPMfilter = data.frame(),
    initialize = function(peaksetfnm, ncell, nfragm) {
      # assume peakSetfnm alwasy has colnames
      self$peakset <- data.table::fread(
        file = peaksetfnm, header = TRUE,
        data.table = FALSE
      )
      self$ncell <- ncell
      self$nfragm <- nfragm
      self$minMacs2Score <- getLog10NegQvalThreshold(nfragm)
    },
    filterByQval = function() {
      message("Soft filtering with nlog10q: ", self$minMacs2Score)
      cond <- self$peakset$score >= self$minMacs2Score
      message(sum(cond), " of ", nrow(self$peakset), " peak left.")
      self$peaksetSoftfilter <- self$peakset[cond, ]
      invisible(self)
    },
    filterBySPM = function(spm = 2) {
      message("SPM filtering with spm: ", spm)
      self$minSPM <- spm
      cond <- self$peakset$spm >= spm
      message(sum(cond), " of ", nrow(self$peakset), " peak left.")
      self$peaksetSPMfilter <- self$peakset[cond, ]
      invisible(self)
    }
  ))

l4s <- L4pcMeta$pL4
peakset_dir <- file.path(work_dir, "out/tscc/itermerge")

#' @param l4 is L4 cluster for peak calling
initL4PeakSet <- function(l4) {
  modifiedL4 <- L4pcMeta[l4, "gL4"]
  peaksetfnm <- file.path(peakset_dir,
    str_glue("{modifiedL4}.filterNfixed.peakset"))
  ncell <- L4pcMeta[l4, "size"]
  nfragm <- L4pcMeta[l4, "rmillion"]
  r <- L4PeakSet$new(peaksetfnm, ncell, nfragm)
  r$filterByQval()
  r$filterBySPM(spm = 2)
  return(r)
}

# This is faster than expected
l4peaksets <- map(.x = l4s, .f = \(l4) {
  modifiedL4 <- L4pcMeta[l4, "gL4"]
  peaksetfnm <- file.path(peakset_dir,
    str_glue("{modifiedL4}.filterNfixed.peakset"))
  ncell <- L4pcMeta[l4, "size"]
  nfragm <- L4pcMeta[l4, "rmillion"]
  r <- L4PeakSet$new(peaksetfnm, ncell, nfragm)
  r$filterByQval()
  r$filterBySPM(spm = 2)
  return(r)
})
names(l4peaksets) <- l4s
# too big
# saveRDS(l4peaksets, "tmp.L4PeakSets.rds")


# * after filtering, merge all the peaks.
# read bed to gr
read2gr <- function(bedF, label) {
  df <- fread(bedF, sep = "\t", header = F)
  colnames(df) <- c("chr", "start", "end", "name", "score")
  df$label <- label
  gr <- GRanges(
    df$chr,
    IRanges(df$start, df$end)
  )
  mcols(gr)$score <- df$score
  mcols(gr)$name <- df$name
  mcols(gr)$label <- df$label
  return(gr)
}

# get non-overlapped regions ' Retreive a non-overlapping set of
# regions from a Genomic Ranges object
#'
#' This function returns a GRanges object containing a non-overlapping
#' set regions derived from a supplied Genomic Ranges object.
#'
#' @param gr A `GRanges` object.
#' @param by The name of a column in `mcols(gr)` that should be used
#'   to determine how overlapping regions should be resolved.  The
#'   resolution of overlapping regions also depends on
#'   `decreasing`. For example, if a column named "score" is used for
#'   `by`, `decreasing = TRUE` means that the highest "score" in the
#'   overlap will be retained and `decreasing = FALSE` means that the
#' lowest "score" in the overlap will be retained.
#' @param decreasing A boolean value indicating whether the values in
#'   the column indicated via `by` should be ordered in decreasing
#'   order. If `TRUE`, the higher value in `by` will be retained.
#' @param verbose
#' A boolean value indicating whether the output should include extra reporting.
#' @export
nonOverlappingGR <- function(
    gr = NULL, 
    by = "score", 
    decreasing = TRUE, 
    verbose = FALSE
  ) {
  stopifnot(by %in% colnames(mcols(gr)))
  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  .clusterGRanges <- function(gr = NULL, filter = TRUE,
                              by = "score", decreasing = TRUE) {
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  if(verbose) {
    message("Converging", appendLF = FALSE)
  }
  i <-  0
  grConverge <- gr
  while (length(grConverge) > 0) {
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i <-  i + 1
    grSelect <- .clusterGRanges(
      gr = grConverge, 
      filter = TRUE, 
      by = by, 
      decreasing = decreasing)

    grConverge <- subsetByOverlaps(
      grConverge,
      grSelect, 
      invert=TRUE, 
      ignore.strand = TRUE) #blacklist selected gr
    if(i == 1){
      #if i=1 then set gr_all to clustered
      grAll <- grSelect
    }else{
      grAll <- c(grAll, grSelect)
    } 
  }
  message(sprintf("Converged after %s iterations!", i))
  if(verbose){
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))
  return(grAll)
}
## merged.gr <- nonOverlappingGR(merged.gr, by = "spm", decreasing = TRUE)
## outUnion <- as.data.frame(merged.gr)


# * check ovlp with our sa1
# * check ovlp with cembav1
# * check ovlp with encoder rDHS
# * homer annotation
# * prepare pmat
# * prepare random pmat
# * add saturation analysis
