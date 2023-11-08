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
sa2L4SumDir <- file.path(here::here(),
  "17.snapatac2/result", "L4_sa2_dlt2_pkl")

files <- list.files(path = sa2L4SumDir, pattern = "\\.pkl$")
L3Clusters <- map_chr(files, function(f) {
  str_split_1(f, "[_.]")[4]
})
fullfiles <- paste(sa2L4SumDir, files, sep = "/")

nm_prefix <- "sa2_clustering"
cll <- "L3"
leidenSums <- purrr:::map(L3Clusters, function(i) {
  message("Loading cluster:", i, "'s LeisenSum.")
  f <- file.path(sa2L4SumDir, str_glue("{nm_prefix}_{cll}_{i}.pkl"))
  newLeidenSumFromPickleFile(pickleFile = f,
    cll = cll, cid = i, nm_prefix = nm_prefix,
    cidIsInt = FALSE,
    cidFromZero = FALSE,
    python = "/home/szu/mambaforge/envs/sa2/bin/python")
})
saveRDS(object = leidenSums, file = "sa2L4LeidenSum.rds")

# * resolution selection per L3 cluster
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

orderstrby3 <- function(arrayOfStr, split = "-") {
  a_b_c <- map(arrayOfStr, function(i) {
    map_int(stringr::str_split_1(i, split), as.integer)
  }) |> do.call(what = rbind)
  order(a_b_c[, 1], a_b_c[, 2], a_b_c[,3])
}

## //tmp comment region below to escape syntax highlight bug
arrayOfStr <- rByMaxSils |>
  _$name |>
  x => gsub(stringr::str_glue("{nm_prefix}_{cll}_"), "", x)
## //

rByMaxSils2 <- orderstrby3(arrayOfStr, split = "-") |>
  x => rByMaxSils[x, ]

write.table(rByMaxSils2, file = "rByMaxSilsL4.csv",
  quote = F, sep = ",", row.names = F, col.names = F)

## for fix order issue of rByMaxSils
## rByMaxSils <- read.table("rByMaxSilsL4.csv", header = F, sep = ",")
## colnames(rByMaxSils) <- c("ith", "sil", "nc", "size", "name")
## add one more column for later manual check label
tmp <- read.table("rByMaxSilsL4.csv", header = F, sep = ",")
tmp$modify <- 0
write.table(tmp, file = "rByMaxSilsL4_v2.csv",
  quote = F, sep = ",", row.names = F, col.names = F)
## for fix order issue of rByMaxSils
rByMaxSils <- read.table("rByMaxSilsL4_v2.csv", header = F, sep = ",")
colnames(rByMaxSils) <- c("ith", "sil", "nc", "size", "name", "modify")
arrayOfStr <- rByMaxSils |>
  _$name |>
  x => gsub(stringr::str_glue("{nm_prefix}_{cll}_"), "", x)
## //
rByMaxSils2 <- orderstrby3(arrayOfStr, split = "-") |>
  x => rByMaxSils[x, ]
write.table(rByMaxSils2, file = "rByMaxSilsL4_v3.csv",
  quote = F, sep = ",", row.names = F, col.names = F)
# we later use rByMaxSilsL4_v3.csv  to replace rByMaxSilsL4.csv


# * add all the clusters
L3Reso <- cluSumBySa2$loadL3Resos()
barcode2L3 <- cluSumBySa2$loadbarcode2L3()
uL3s <- unique(barcode2L3$L3)
L3inL4 <- gsub("sa2_clustering_L3_", "", rByMaxSils$name)
totclu <- sum(rByMaxSils2$nc) + (length(uL3s)-length(L3inL4))



# manually check
# load manually checked file and then double check annots
L3Sums <- cluSumBySa2$loadL3Sums()
L3Resos <- cluSumBySa2$loadL3Resos()
barcode2L3 <- cluSumBySa2$loadbarcode2L3()

L4Resos <- cluSumBySa2$loadL4Resos()
L4Sums <- cluSumBySa2$loadL4Sums()


tmp <- map_lgl(L4Sums, function(i) {
  cid <- i$cid
  rowOfReso <- L4Resos[paste0("sa2v1_", cid), ]
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

## # all pass

## # * summarize barcode to clustering results
barcode2L4List <- map(L4Sums, function(i){
  cid <- i$cid
  message("Current cid: ", cid)
  rowOfReso <- L4Resos[paste0("sa2v1_", cid), ]
  index <- rowOfReso$whichMaxSil
  r <- i$rs[index]
  i$getBarcode2ClusterName(r = r, nmSep = "-", cll = "L4")
})

barcode2L4 <- do.call(rbind, barcode2L4List)

# append L3s' that are not in L4s.
barcode2L3_left <- subset(barcode2L3, !(barcode %in% barcode2L4$barcode))
barcode2L3_left$L4 <- paste(barcode2L3_left$L3, "1", sep = "-")
barcode2L3_left$L3 <- NULL

# concat barcode2L4 and barcode2L3_left to barcode2L4_final
barcode2L4_final <- rbind(barcode2L4, barcode2L3_left)
rownames(barcode2L4_final) <- barcode2L4_final$barcode

# organize barcode2L4_final based on L4
# use vapply is faster than original funciton map and do.call
asIntL4 <- vapply(barcode2L4_final$L4, function(i) {
  stringr::str_split(i, "-")[[1]]
}, rep("1", 4)) |> as.data.frame(t())
asIntL4_bk <- asIntL4
asIntL4 <- t(asIntL4)
# transform to integer
asIntL4 <- asIntL4 |>
  x => apply(x, 2, as.integer) |>
  as.data.frame()
# sort by column
index <- order(asIntL4[,1], asIntL4[,2], asIntL4[,3], asIntL4[,4])
barcode2L4_sort <- barcode2L4_final[index,]
rownames(barcode2L4_sort) <- barcode2L4_sort$barcode

# check if barcode2L4_sort is same as barcode2L4_final
a1 <- table(barcode2L4_sort$L4)
a2 <- table(barcode2L4_final$L4)
sum(a1 == a2)

sum(barcode2L4_sort[barcode2L4$barcode, "L4"] == barcode2L4$L4)
sum(barcode2L4_sort[barcode2L3_left$barcode, "L4"] == barcode2L3_left$L4)

write.table(barcode2L4_sort,
  file = "sa2.barcode2L4.csv", quote= FALSE, sep = ",",
  col.names = TRUE, row.names = FALSE)

## # check barcode2L4 file
barcode2L4 <- cluSumBySa2$loadbarcode2L4()
## # * check overall clustering stats
clusters <- str_split(barcode2L4$L4, "-", simplify = TRUE)
table(clusters[, 1])
length(table(paste(clusters[,1], clusters[,2], sep = "-")))
length(table(barcode2L4$L4))
# * select some cluster id, then check barcodes are matched.
L4Sums <- cluSumBySa2$loadL4Sums()
L4Resos <- cluSumBySa2$loadL4Resos()
cid <- "10-1-1"
l <- L4Sums[[cid]]
s <- L4Resos[str_glue("sa2v1_{cid}"), ]
barcode2L4_l <- data.frame(
  barcode = l$barcodes,
  L4  = paste(cid, l$label[, s[1,"whichMaxSil"]], sep = "-")
)
all(barcode2L4_l$L4 == barcode2L4[barcode2L4_l$barcode, "L4"])
