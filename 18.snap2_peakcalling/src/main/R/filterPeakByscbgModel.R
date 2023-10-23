library(tidyverse)
library(rlang)
library(rlang)
library(data.table)
library(reticulate)
library(furrr)
plan(multicore, workers = 5)
options(future.globals.maxSize= 15e9)
library(future.apply)

library(gamlss.dist)
library(gamlss)
library(fitdistrplus)

rdir <- file.path(here::here(), "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env, Sa2PeakCalling)

usePython <- "/home/szu/mambaforge/envs/sa2/bin/python"
reticulate::use_python(usePython)
pd <- reticulate::import("pandas", convert = FALSE)

# * meta
atacMeta <- readRDS(cembav2env$sa2metaFile)
projdir <- here::here()
workdir <- file.path(projdir, "18.snap2_peakcalling")
default_cutoff <- 0.001
pep <- 0.001
rnd_upbound <- 0.1
qval <- 0.01
outdir <- file.path(workdir, "out/scfilter")

fitbgdir <- file.path(workdir, "out/scfilter",
  "fitfrac_bg")
fitbgsuffix <- ".fitPeakModel.para.csv"

loadfitbg <- function(l4) {
  data.table::fread(
    file = file.path(fitbgdir, paste0(l4, fitbgsuffix)),
    sep = ",",
    header = TRUE,
    data.table = FALSE)
}


# * read union peak
peakfrac_union <- pd$read_pickle(
  file.path(workdir, "out/scfilter", "peakfrac_union.pkl"))
peakfrac <- reticulate::py_to_r(x = peakfrac_union)

# * check all the models are converged
# all of them converged
pL4s <- rownames(peakfrac)
checkConverge <- future_map_lgl(pL4s,
  \(x) {
    m <- loadfitbg(pL4)
    if(m$converged > 0) {
      TRUE
    } else {
      FALSE
    }
  }
)

# * fit model per pL4
pvalues <- future_map(.x = pL4s, .f = \(l4) {
  m <- loadfitbg(l4)
  p <- 1 - gamlss.dist::pBEZI(
    unlist(peakfrac[l4, ]), mu = m$mu, sigma = m$sigma,
    nu = m$nu, lower.tail = TRUE, log.p = FALSE)
  return(p)
})
names(pvalues) <- pL4s

# check pvalues
## (pL4 <- sample(pL4s, 1))
## p1 <- pvalues[[pL4]]
## m <- loadfitbg(pL4)
## t1 <- 1-gamlss.dist::pBEZI(
##   unlist(peakfrac[pL4, ]), mu = m$mu, sigma = m$sigma,
##   nu = m$nu, lower.tail = TRUE, log.p = FALSE)
## sum(p1 == t1) == length(p1)

## save pvalues as matrix
pvaluesMat <- do.call(rbind, pvalues)
rownames(pvaluesMat) <- pL4s
colnames(pvaluesMat) <- colnames(peakfrac)
# > 6G
saveRDS(pvaluesMat,
  file = file.path(outdir, "pL4.unionpeak.pvalue.mat.rds"))

# * qvalue and filter
unionpeaks <- Sa2PeakCalling$readUnionPeakSetFile()
spm <- 5
filterUnionPeaks <- rownames(unionpeaks)[unionpeaks$spm >= spm]
filterPvalMat <- pvaluesMat[, filterUnionPeaks]

qvaluesUnion <- vapply(pL4s, FUN = \(l4) {
  message(l4)
  p.adjust(filterPvalMat[l4, ], method = "BH")
}, FUN.VALUE = filterPvalMat[1, ])

qvaluesUnion <- t(qvaluesUnion)

min_qval <- furrr::future_map_dbl(.x = seq_len(ncol(qvaluesUnion)),
  .f = \(i) {
    min(qvaluesUnion[, i])
  }
)
names(min_qval) <- colnames(qvaluesUnion)

scfilteredPeaks <- names(min_qval)[min_qval <= qval]
## all the peaks pass the sc-filtering

# check all union peak
qvaluesAll <- future.apply::future_vapply(
  pL4s, FUN = \(l4) {
    message(l4)
    p.adjust(pvaluesMat[l4, ], method = "BH")
  }, FUN.VALUE = pvaluesMat[1, ])
colnames(qvaluesAll) <- pL4s

## all the peaks pass the sc-filtering
min_qval_all <- furrr::future_map_dbl(.x = seq(nrow(qvaluesAll)),
  .f = \(i) {min(qvaluesAll[i, ])})
sum(min_qval_all <= qval)

# * test
## i <- 1
## pL4 <- rownames(peakfrac)[i]
## m <- loadfitbg(pL4)
## p <- 1 - gamlss.dist::pBEZI(
##   unlist(peakfrac[i, ]), mu = m$mu, sigma = m$sigma,
##   nu = m$nu, lower.tail = TRUE, log.p = FALSE)
## q <- p.adjust(p, method = "BH")

