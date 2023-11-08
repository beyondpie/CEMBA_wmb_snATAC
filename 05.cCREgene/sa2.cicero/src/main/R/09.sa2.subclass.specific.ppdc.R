library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

projdir <- here::here()
rpacdir <- file.path(projdir, "package/R")
workdir <- file.path(projdir, "04.cCREgene/sa2.cicero")
pdcsumdir <- file.path(workdir, "out/sa2pdcsum")
pdcdir <- file.path(workdir, "out/sa2pdc")

# * meta
subclassfnm <- file.path(projdir, "meta",
  "sa2.subclass.src.txt")

# * reduce redundant in pos.pdc.bedpe
# * functions
read_ppdc_bedpe10 <- function(fnm) {
  nms <- c("pchr", "pstart", "pend", "dchr", "dstart", "dend",
    "pair", "coacc", "n1", "n2")
  data.table::fread(file = fnm, sep = "\t", header = FALSE,
    data.table = FALSE) |>
    setNames(object = _, nm = nms)
}

uniq_ppdc_bedpe_fn <- function(fnm, outfnm) {
  rawppdc_bedpe <- read_ppdc_bedpe10(fnm)
  uniq_ppdc_bedpe <- rawppdc_bedpe |>
    group_by(pair) |>
    slice_max(order_by = coacc, n = 1)
  write.table(uniq_ppdc_bedpe, outfnm,
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = FALSE)
  return(uniq_ppdc_bedpe)
}

# * uniqe pdc bedpe
## pos pdc
uniq_ppdc_bedpe_fn(
  fnm = file.path(pdcsumdir,
  "mba.whole.sa2subclass.pearson.pos.pdc.bedpe"),
  outfnm = file.path(workdir, "out", "sa2pdcsum",
    "mba.whole.sa2subclass.pos.pdc.uniq.coacc.bedpe"))

## neg pdc
uniq_ppdc_bedpe_fn(
  fnm = file.path(pdcsumdir,
  "mba.whole.sa2subclass.pearson.neg.pdc.bedpe"),
  outfnm = file.path(workdir, "out", "sa2pdcsum",
    "mba.whole.sa2subclass.neg.pdc.uniq.coacc.bedpe"))

# uniq rdm
uniq_ppdc_bedpe_fn(
  fnm = file.path(pdcsumdir,
  "sa2subclass.pdc.rdm.dup.bedpe"),
  outfnm = file.path(workdir, "out", "sa2pdcsum",
    "mba.whole.sa2subclass.rdm.pdc.uniq.coacc.bedpe"))

# * get subclass-specific ppdc
out_subclass_ppdc_dir <- file.path(workdir, "out/sa2ppdc")
if (!dir.exists(out_subclass_ppdc_dir)) {
  dir.create(out_subclass_ppdc_dir)
}


all.ppdc <- read_ppdc_bedpe10(
  fnm = file.path(pdcsumdir,
    "mba.whole.sa2subclass.pos.pdc.uniq.coacc.bedpe")
)

asct.pdc <- read_ppdc_bedpe10(
  fnm = file.path(pdcdir, "sa2subclass.Astro-TE_NN.pdc.bedpe")
)

asct.ppdc <- asct.pdc[asct.pdc$pair %in% all.ppdc$pair, ] |>
  group_by(pair) |>
  slice_max(order_by = coacc, n = 1)
write.table(asct.ppdc, file.path(out_subclass_ppdc_dir,
  "sa2subclass.Astro-TE_NN.pos.pdc.bedpe"), quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = FALSE)


subs <- data.table::fread(
  file = file.path(here::here(), "23.m3C", "subclass.txt"),
  header = FALSE, data.table = FALSE)$V1

for (i in seq_along(subs)) {
  sc <- subs[i]
  message("processing ", sc)
  scpdc <- read_ppdc_bedpe10(
    fnm = file.path(pdcdir, paste0("sa2subclass.", sc, ".pdc.bedpe"))
  )
  scpdc <- scpdc[scpdc$pair %in% all.ppdc$pair, ] |>
    group_by(pair) |>
    slice_max(coacc, n = 1)

  scpdc <- scpdc[order(scpdc$coacc, decreasing = TRUE), ]
  top20_scpdc <- scpdc[1:round(nrow(scpdc) * 0.20), ]
  end20_scpdc <- scpdc[(nrow(scpdc) - round(nrow(scpdc) * 0.80)):nrow(scpdc), ]

  write.table(top20_scpdc, file.path(out_subclass_ppdc_dir,
    paste0("sa2subclass.", sc, ".top20.ppdc.bedpe")), quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = FALSE)

  write.table(end20_scpdc, file.path(out_subclass_ppdc_dir,
    paste0("sa2subclass.", sc, ".end20.ppdc.bedpe")), quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = FALSE)
}
