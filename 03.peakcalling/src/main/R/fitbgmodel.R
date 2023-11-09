library(tidyverse)
library(reticulate)
library(gamlss.dist)
library(gamlss)
library(fitdistrplus)

usePython <- "/home/szu/mambaforge/envs/sa2/bin/python"
reticulate::use_python(usePython)
pd <- reticulate::import("pandas", convert = FALSE)

projdir <- here::here()
workdir <- file.path(projdir, "18.snap2_peakcalling")

# * meta
default_cutoff <- 0.001
pep <- 0.001
rnd_upbound <- 0.1
outdir <- file.path(workdir, "out/scfilter/fitfrac_bg")


# * read peakfrac data
peakfrac_rnd <- pd$read_pickle(
  file.path(workdir, "out/scfilter", "peakfrac_rnd.pkl"))
peakfrac_rnd <- reticulate::py_to_r(x = peakfrac_rnd)

## peakfrac_union <- pd$read_pickle(
##   file.path(workdir, "out/scfilter", "peakfrac_union.pkl"))
## peakfrac_union <- reticulate::py_to_r(x = peakfrac_union)

# * fit models for each pL4 using peakfrac_rnd
i <- as.integer(commandArgs(trailingOnly = TRUE)[1])
pL4 <- rownames(peakfrac_rnd)[i]
message("fit BEZI model for ", pL4)

x <- unlist(peakfrac_rnd[i, ])
x <- x[x <= rnd_upbound]

mod <- gamlss::gamlss(
  x ~ 1, sigma.formula = ~1,
  nu.formula = ~1,
  family = BEZI,
  control = gamlss::gamlss.control(n.cyc = 100, trace = FALSE),
  )

## summary(mod)
mufit <- fitted(mod, "mu")[1]
sigmafit <- fitted(mod, "sigma")[1]
nufit <- fitted(mod, "nu")[1]
converged <- ifelse(mod$converged, 1, 0)

cutoff <- if (mod$converged) {
  # find x axis value at p = pep (default 0.001)
  gamlss.dist::qBEZI((1 - pep),
    mu = mufit, sigma = sigmafit, nu = nufit,
    lower.tail = TRUE, log.p = FALSE)
} else {
  message("Fitting BEZI dose not converge.",
    " Use default cutoff: ", default_cutoff)
  default_cutoff
}

para <- data.frame(
  n = length(x),
  mu = mufit,
  sigma = as.numeric(sigmafit),
  nu = nufit,
  pep = pep,
  cutoff = cutoff,
  Gdeviance = mod$G.deviance,
  converged = converged)

data.table::fwrite(para,
  file.path(outdir, paste0(pL4, ".fitPeakModel.para.csv")),
  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")

message("done.")
