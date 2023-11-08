library(data.table)
library(purrr)
library(Matrix)

packdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env, CEMBAMeta)

# TODO: put this class to R package
L3ClusterPeakSum <- R6::R6Class(
  classname = "L3ClusterPeakSum",
  public = list(
    L3Id = "L3id",
    name = "L3Annot",
    subclass = "allen_annot",
    # index of 1.3 millio peaks we have
    peaks = c(1, 2),
    regions = c(0.1),
    barcodes = c("CEMBA1", "CEMBA2"),
    initialize = function(L3Id,
                           name,
                           # binary pmat
                           bpmat,
                           # barcode to L3
                           b2L3df,
                           regions,
                           # barcode 2 region
                           b2rdf,
                           # L3 to subclass
                           L3toSubclass) {
      self$L3Id <- L3Id
      self$name <- name
      self$subclass <- L3toSubclass[self$L3Id, 2]
      self$peaks <- which(bpmat[self$L3Id, ] == 1)
      self$barcodes <- b2L3df[b2L3df[, 2] == self$L3Id, 1]
      self$regions <- getRegions(b2rdf, self$barcodes, regions)
    },
    getMainRegion = function(threshold = 0.5) {
      names(self$regions)[self$regions >= threshold]
    },
    isMainRegion = function(region = "HY", threshold = 0.5) {
      mr <- self$getMainRegion(threshold)
      ifelse(region %in% mr, TRUE, FALSE)
    }
    # TODO: add print function
  )
)

# TODO: merge to cembaObject
getRegions <- function(b2rdf, barcodes, regions) {
  r <- seq_len(length(regions))
  names(r) <- regions
  cnts <- table(b2rdf[barcodes, 2])
  r[names(cnts)] <- cnts
  r[! (names(r) %in% names(cnts))] <- 0
  r / sum(r)
}

# * load data
cembaMeta <- readRDS(cembav2env$allmetaFile)
cembaMeta$barcodes <- with(cembaMeta, paste(sample, barcode, sep = "."))

cembaObject <- CEMBAMeta$new()
allPeakInfo <- cembaObject$readAllPeaks(cembav2env)
bpmatL3 <- cembaObject$readL3BinaryPmat(cembav2env)
majorRegions <- unique(cembaMeta$MajorRegion)

# TODO: merge these info to cembaObject
L3s <- rownames(bpmatL3)
L3toL3Annot <- unique(cembaMeta[, c("L3", "upL3Annot")])
rownames(L3toL3Annot) <- L3toL3Annot$L3
L3toSubclass <- unique(cembaMeta[, c("L3", "AllenAnnotConcat")])
rownames(L3toSubclass) <- L3toSubclass$L3

b2L3df <- cembaMeta[, c("barcodes", "L3")]
b2rdf <- cembaMeta[, c("barcodes", "MajorRegion")]

# * get HY-majored peaks
mr <- "HY"
l3PeakSumList <- purrr::map(L3s, function(i) {
  print("generate peak sum for: ", i)
  L3ClusterPeakSum$new(
    L3Id = i,
    name = L3toL3Annot[i, 2],
    bpmat = bpmatL3,
    b2L3df = b2L3df,
    regions = majorRegions,
    b2rdf = b2rdf,
    L3toSubclass = L3toSubclass
  )
})
names(l3PeakSumList) <- L3s
L3s_HY <- L3s[purrr::map_lgl(L3s, function(i) {
  l3PeakSumList[[i]]$isMainRegion("HY", threshold = 0.55)
})]
L3PeakSum.HYMajor <- l3PeakSumList[L3s_HY]
saveRDS(L3PeakSum.HYMajor, "L3PeakSum.HYMajor.rds")
allHYPeakIds <- unique(unlist(purrr::map(
  L3PeakSum.HYMajor, function(i) {i$peaks}
))) |> sort()
allHYPeaks <- allPeakInfo$peaks[allHYPeakIds, ]
write.table(allHYPeaks[, c("chr", "start", "end")],
  file = "peaks.HY.associated.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE)

# * test
barcodes <- with(cembaMeta, barcodes[L3 == "GABA_4_9"])
r <- getRegions(b2rdf = b2rdf,
  barcodes = barcodes, regions = majorRegions)
a <- L3ClusterPeakSum$new(
  L3Id = "GABA_4_9",
  name = L3toL3Annot["GABA_4_9", 2],
  bpmat = bpmatL3,
  b2L3df = b2L3df,
  regions = majorRegions,
  b2rdf = b2rdf,
  L3toSubclass = L3toSubclass
)
