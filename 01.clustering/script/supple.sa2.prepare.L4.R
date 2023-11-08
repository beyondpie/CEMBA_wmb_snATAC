library(data.table)
library(stringr)
library(dplyr)
packdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cluSumBySa2)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

# * load L3 clustering result
L3Sums <- cluSumBySa2$loadL3Sums()
L3Resos <- cluSumBySa2$loadL3Resos()
barcode2L3 <- cluSumBySa2$loadbarcode2L3()

# * prepare two files
# 0. get L3 needed for L4 clustering
uL3s <- unique(barcode2L3$L3)

# nolint start
L3Pattern <- with(L3Resos, L2nm[needL4 == "L4"]) |>
  x => gsub("sa2v1_", "", x)
nMin <- 400
# nolint end

# filter uL3s that are in L3Pattern
uL3toL2 <- gsub("-[0-9]+$", "", uL3s)
uL3forL4 <- uL3s[uL3toL2 %in% L3Pattern]

# 1. cluster size at L3-level in uL3forL4
# and no less than nMin
L3tosize <- table(barcode2L3$L3) |> 
  as.data.frame(stringsAsFactors = FALSE) |>
  setNames(c("L3", "size")) |>
  x => x[x$L3 %in% uL3forL4, ] |>
  arrange(desc(size)) |>
  x => x[x$size >= nMin, ]
# outut L3tosize to csv file without header
write.table(L3tosize,
  file = file.path("../resource", "sa2_dlt2_L3forL4_cluster2size.csv"),
  sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 2. barcode2cluster at L3-level
barcode2L3 |> 
  x => x[x$L3 %in% uL3forL4, ] |> 
  write.table(
    file = file.path("../resource", "sa2_dlt2_L3forL4_barcode2id.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

