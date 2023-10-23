# Use the rough annotation from snATACv1 to annotate the snATACv2 data.
# especially check IMN, NN, and CB granule cells.
# this is especially commented in snmc manuscript.
library(dplyr)
library(stringr)
library(tibble)
library(purrr)
library(furrr)
library(progressr)
packdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env, cembaMeta.sa1, cluSumBySa2)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
plan(multicore, workers = 5)
# 8G for future
options(future.globals.maxSize = 8000 * 1024^2)
# * load snATAc meta of SnapATAC1
barcode2annot.sa1 <- cembaMeta.sa1$getBarcode2Annot.sa1(cembav2env)
barcode2clu.sa2 <- cluSumBySa2$loadbarcode2L3()

# * generate rough annotations based on snATACv1
qualityMap <- data.frame(
  quality = c("low", "medium", "high"),
  minFrac = c(0.0, 0.4, 0.7)
)

barcode.joint <- intersect(barcode2annot.sa1$barcode,
  barcode2clu.sa2$barcode)
barcode2annot.joint <- barcode2annot.sa1[barcode.joint, ]

uL3s.sa2 <- unique(barcode2clu.sa2$L3)

mapCluBySa2ToAnnotSa1 <- function(clu.sa2,
                                  annot.level = "AllenAnnot", top = 1) {
  rawBarcodes <- barcode2clu.sa2[barcode2clu.sa2[,2] == clu.sa2, 1]
  barcodes <- intersect(rawBarcodes, barcode.joint)
  if (length(barcodes) == 0) {
    message("no barcodes detected for ", clu.sa2, " in snATAC.sa1.")
    return(c(clu.sa2, "0", "0", "0", NA, NA))
  }
  r <- length(barcodes) / length(rawBarcodes)
  s <- table(barcode2annot.sa1[barcodes, annot.level])
  sRatio <- s/sum(s)
  s <- sRatio |>
    as.data.frame(stringsAsFactors = FALSE) |>
    setNames(c(annot.level, "ratio")) |>
    dplyr::arrange(dplyr::desc(ratio)) |>
    head(top)
  quality <- if (s$ratio < qualityMap$minFrac[2]) {
    qualityMap$quality[1]
  } else if (s$ratio < qualityMap$minFrac[3]) {
    qualityMap$quality[2]
  } else {
    qualityMap$quality[3]
  }
  return(c(
    clu.sa2,
    as.character(length(barcodes)),
    format(r, format = "f", digits = 2),
    s[[annot.level]],
    quality,
    format(s$ratio, format = "f", digits = 2)
    )
  )
}

mapAllClusToAnnotSa1.async <- function(annot.level = "AllenAnnot",
                                       clus = uL3s.sa2, top = 1) {
  with_progress({
    p <- progressor(steps = length(clus))
    r <- future_map(clus, \(x) {
      p()
      mapCluBySa2ToAnnotSa1(x, annot.level, top = 1)
    }) |>
      # nolint start
      x => do.call(rbind, x) |>
      # nolint end
    as.data.frame(stringsAsFactors = FALSE) |>
    setNames(
      c("clu.sa2",
      "size",
      "sa1_ratio",
      annot.level,
      paste0(annot.level, "_quality"),
      paste0(annot.level, "_ratio")
      ))
  })
  rownames(r) <- r$clu.sa2
  return(r)
}

L3toAllenAnnot.rough.sa2 <- mapAllClusToAnnotSa1.async("AllenAnnot")
L3toL3Annot.rough.sa2 <- mapAllClusToAnnotSa1.async("L3Annot")
L3toL2Annot.rough.sa2 <- mapAllClusToAnnotSa1.async("L2Annot")
L3toclass.rough.sa2 <- mapAllClusToAnnotSa1.async("class")
L3tomainclass.rough.sa2 <- mapAllClusToAnnotSa1.async("mainclass")

# merge about data.frame with column: clu.sa2
L3toSa1Annot.rough.sa2 <- list(L3toAllenAnnot.rough.sa2,
  L3toL3Annot.rough.sa2,
  L3toL2Annot.rough.sa2,
  L3toclass.rough.sa2,
  L3tomainclass.rough.sa2) |>
  reduce(left_join, by = c("clu.sa2", "size", "sa1_ratio"))
rownames(L3toSa1Annot.rough.sa2) <- L3toSa1Annot.rough.sa2$clu.sa2

# save to csv
write.table(L3toSa1Annot.rough.sa2,
  file = file.path("../resource", "sa2_dlt2_L3toSa1Annot.rough.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# * Annotate snATACv2 data
barcode2clu2annot <- barcode2clu.sa2
barcode2clu2annot$mainclass <- L3toSa1Annot.rough.sa2[barcode2clu.sa2$L3, "mainclass"]
barcode2clu2annot$class <- L3toSa1Annot.rough.sa2[barcode2clu.sa2$L3, "class"]
barcode2clu2annot$AllenAnnot <- L3toSa1Annot.rough.sa2[barcode2clu.sa2$L3, "AllenAnnot"]

# * save to csv
write.table(barcode2clu2annot,
  file = file.path("../resource", "sa2_dlt2_barcode2L3clu2annot.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
