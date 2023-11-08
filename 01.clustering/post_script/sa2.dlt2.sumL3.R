library(stringr)
library(ggplot2)
library(purrr)
library(dplyr)

packdir <- file.path(here::here(), "package/R")
import::from(.from = "sa2Leiden.R", .directory = packdir,
  LeidenSum,
  newLeidenSumFromPickleFile,
  getCluster2Size,
  getBarcode2ClusterName)

import::from(.from = "cembav2env.R", .directory = packdir,
  cluSumBySa2)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * load L3 LeidenSum
sa2L3SumDir <- file.path(here::here(),
  "17.snapatac2/result", "L3_sa2_dlt2_pkl")
files <- list.files(path = sa2L3SumDir, pattern = "\\.pkl$")
L2Clusters <- map_chr(files, function(f) {
  str_split_1(f, "[_.]")[4]
})
fullfiles <- paste(sa2L3SumDir, files, sep = "/")

nm_prefix <- "sa2_clustering"
cll <- "L2"
leidenSums <- purrr:::map(L2Clusters, function(i) {
  message("Loading cluster:", i, "'s LeisenSum.")
  f <- file.path(sa2L3SumDir, str_glue("{nm_prefix}_{cll}_{i}.pkl"))
  newLeidenSumFromPickleFile(pickleFile = f,
    cll = cll, cid = i, nm_prefix = nm_prefix,
    cidIsInt = FALSE,
    cidFromZero = FALSE,
    python = "~/mambaforge/bin/python")
})
saveRDS(object = leidenSums, file = "sa2L3LeidenSum.rds")

# * resolution selection per L2 cluster
# by max silhouette
getRByMaxSil <- function(ldSum) {
  ith <- which.max(ldSum$sils)
  sil <- ldSum$rils[ith]
  r <- ldSum$rs[ith]
  nc <- length(unique(ldSum$label[, ith]))
  return(list(ith = ith, sil = sil,
    r = r, nc = nc, size = ldSum$nsample, name = ldSum$name))
}
rByMaxSils <- purrr::map_dfr(leidenSums, getRByMaxSil) |>
  as.data.frame()
rownames(rByMaxSils) <- rByMaxSils$name

orderstrby2 <- function(arrayOfStr, split = "-") {
  a_b <- map(arrayOfStr, function(i) {
    map_int(stringr::str_split_1(i, split), as.integer)
  }) |> do.call(what = rbind)
  order(a_b[, 1], a_b[, 2])
}

## //tmp comment region below to escape syntax highlight bug
## arrayOfStr <- rByMaxSils |>
##   _$name |>
##   x => gsub(stringr::str_glue("{nm_prefix}_{cll}_"), "", x)

## rByMaxSils2 <- orderstrby2(arrayOfStr, split = "-") |>
##   x => rByMaxSils[x, ]

## write.table(rByMaxSils2, file = "rByMaxSils.csv",
##   quote = F, sep = ",", row.names = F, col.names = F)
## //

# manually check
## load manually checked file and then double check annots
L3Sums <- cluSumBySa2$loadL3Sums()
L3Resos <- cluSumBySa2$loadL3Resos()


tmp <- map_lgl(L3Sums, function(i) {
  cid <- i$cid
  rowOfReso <- L3Resos[paste0("sa2v1_", cid), ]
  index <- rowOfReso$whichMaxSil
  if (rowOfReso$reso != i$rs[index]) {
    warning(cid, " reso of ", index, " does not match.")
    return(FALSE)
  }
  if (rowOfReso$ncell != i$nsample) {
    warning(cid, " ncell of ", index, " does not match.")
    return(FALSE)
  }
  ncluFromSum <- length(table(i$label[, index]))
  if (rowOfReso$nclu != ncluFromSum) {
    warning(cid, " nclu of ", index, " does not match.")
    return(FALSE)
  }
  return(TRUE)
})

# all pass

# * summarize barcode to clustering results
barcode2L3List <- map(L3Sums, function(i){
  cid <- i$cid
  message("Current cid: ", cid)
  rowOfReso <- L3Resos[paste0("sa2v1_", cid), ]
  index <- rowOfReso$whichMaxSil
  r <- i$rs[index]
  i$getBarcode2ClusterName(r = r, nmSep = "-", cll = "L3")
})

barcode2L3 <- do.call(rbind, barcode2L3List)
write.table(barcode2L3,
  file = "sa2.barcode2L3.csv", quote= FALSE, sep = ",",
  col.names = TRUE, row.names = FALSE)

# check barcode2L3 file
barcode2L3 <- cluSumBySa2$loadbarcode2L3()
# * check overall clustering stats
clusters <- str_split(barcode2L3$L3, "-", simplify = TRUE)
table(clusters[, 1])
length(table(paste(clusters[,1], clusters[,2], sep = "-")))
length(table(barcode2L3$L3))
# * select some cluster id, then check barcodes are matched.
L3Sums <- cluSumBySa2$loadL3Sums()
L3Resos <- cluSumBySa2$loadL3Resos()
cid <- "10-1"
l <- L3Sums[[cid]]
s <- L3Resos[str_glue("sa2v1_{cid}"), ]
barcode2L3_l <- data.frame(
  barcode = l$barcodes,
  L3  = paste(cid, l$label[, s[1,"whichMaxSil"]], sep = "-")
)
all(barcode2L3_l$L3 == barcode2L3[barcode2L3_l$barcode, "L3"])
