library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(future)
library(furrr)
plan(multicore, workers = 7)
options(future.globals.maxSize = 15e9)

library(GenomicRanges)
library(rtracklayer)

projdir <- here::here()
workdir <- file.path(projdir, "04.cCREgene/sa2.cicero")
packdir <- file.path(here::here(), "package", "R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env, Sa2PeakCalling)
import::from(.from = "bed.R", .directory = packdir,
  loadGRFromBed, getGeneTSS, resizeGeneTSS, getGeneGR)

# * configs
window.size <- 2*500*1000
peakfile <- Sa2PeakCalling$finalpeakBedFile
allpdcPairfnm <- cembav2env$sa2.all.pdc.pairs.only
# * load data
allPeaks <- data.table::fread(
  file = peakfile, sep = "\t", header = FALSE,
  data.table = FALSE)
allPeakGR <- loadGRFromBed(bedFile = peakfile)

all.pdc.pairs <- data.table::fread(
  file = allpdcPairfnm, sep = "\t", header = FALSE,
  data.table = FALSE)
colnames(all.pdc.pairs) <- c("pairs", "proximal")

## mouse gene annotate
gtf <- rtracklayer::import(cembav2env$mm.gencode.vM23.gtf)
genes <- getGeneGR(gtf)
gene.loci <- resizeGeneTSS(gtf, window.size = window.size)

# * for each gene, get peaks without a range around it.
genes.pdc <- unique(vapply(unique(all.pdc.pairs$pairs),
  function(i) {
    unlist(strsplit(i, split = "|", fixed = TRUE))[1]
  }, FUN.VALUE = "a"))
gene.loci.pdc <- gene.loci[match(genes.pdc, mcols(gene.loci)$gene_name)]

distal.peaks.pdc <- unique(vapply(unique(all.pdc.pairs$pairs),
  function(i) {
    unlist(strsplit(i, split = "|", fixed = TRUE))[2]
  }, FUN.VALUE = "peak1"))
distalPeakGR <- allPeakGR[match(distal.peaks.pdc, mcols(allPeakGR)$name)]

ovlp.dg <- findOverlaps(query = distalPeakGR, subject = gene.loci.pdc)

# * sample subset of ramdom gene-peak pairs defined above.
n.dist.on.genes <- table(subjectHits(ovlp.dg))
rdm.dist.on.genes <- furrr::future_map(names(n.dist.on.genes), function(i) {
  distals <- unique(queryHits(ovlp.dg)[subjectHits(ovlp.dg) == as.integer(i)])
  set.seed(2023)
  sample(
    setdiff(seq_len(length(distal.peaks.pdc)), distals),
    length(distals))
})

rdm.dists <- unlist(rdm.dist.on.genes)
rdm.genes <- unlist(lapply(names(n.dist.on.genes), function(i) {
  rep(as.integer(i), n.dist.on.genes[i])
}))

rdm.pair.index <- data.frame(
  subject = rdm.genes,
  query = rdm.dists)
rdm.pair.names <- data.frame(
  gene = mcols(gene.loci.pdc)$gene_name[rdm.genes],
  peak = mcols(distalPeakGR)$name[rdm.dists]
)
saveRDS(rdm.pair.names,
  file = file.path(workdir, "out/sa2pdcsum",
    "sa2.rdm.g2p.pdc.rds"))

# * organize atac pairs
pdc.pair.info <- data.table::fread(
  file = allpdcPairfnm,
  header = FALSE, data.table = FALSE, sep = "\t"
)
pdc.pair <- unique(pdc.pair.info$V1)
pdc.pair.df <- t(vapply(pdc.pair, function(i) {
    unlist(strsplit(i, split = "|", fixed = TRUE))
}, FUN.VALUE = c("a", "peak1")))
colnames(pdc.pair.df) <- c("gene", "peak")
rownames(pdc.pair.df) <- NULL
pdc.pair.df <- as.data.frame(pdc.pair.df)
saveRDS(pdc.pair.df,
  file = file.path(workdir, "out/sa2pdcsum",
    "sa2.g2p.pdc.pair.rds"))


# *  get rdm bedpe for hic verification
read_ppdc_bedpe10 <- function(fnm) {
  nms <- c("pchr", "pstart", "pend", "dchr", "dstart", "dend",
    "pair", "coacc", "n1", "n2")
  data.table::fread(file = fnm, sep = "\t", header = FALSE,
    data.table = FALSE) |>
    setNames(object = _, nm = nms)
}
rdm.pdc <- readRDS(file.path(workdir, "out/sa2pdcsum",
  "sa2.rdm.g2p.pdc.rds"))

all.bedpe <- read_ppdc_bedpe10(fnm = file.path(
  here::here(), "04.cCREgene", "sa2.cicero", "out",
  "sa2pdcsum", "mba.whole.sa2subclass.merge.bedpe.all"
))

pdc.pair.df <- t(vapply(all.bedpe$pair, function(i) {
    unlist(strsplit(i, split = "|", fixed = TRUE))
}, FUN.VALUE = c("a", "peak1")))
colnames(pdc.pair.df) <- c("gene", "peak")

g2tss <- data.frame(
  gene = pdc.pair.df[,1],
  chr = all.bedpe$pchr,
  start = all.bedpe$pstart,
  end = all.bedpe$pend
) |> group_by(gene) |> slice_sample(n = 1)

rownames(g2tss) <- g2tss$gene

rdm.pdc <- rdm.pdc[rdm.pdc$gene %in% rownames(g2tss), ]

peakstr2bed <- function(peak) {
  chr2pos <- str_split_1(peak, ":")
  chr <- chr2pos[1]
  pos <- str_split_1(chr2pos[2],"-")
  c(chr,  pos)
}

peak.rmd <- furrr::future_map(rdm.pdc$peak,peakstr2bed)
peak.rdm <- do.call(rbind, peak.rmd)

rdm.p.bed <- g2tss[rdm.pdc$gene, c("chr", "start", "end")]

rdm.bedpe.df <- cbind(rdm.p.bed, peak.rdm)

data.table::fwrite(rdm.bedpe.df,
  file = file.path(here::here(), "04.cCREgene", "sa2.cicero",
    "out", "sa2pdcsum", "sa2subclass.rdm.g2p.bedpe"),
  row.names = FALSE, col.names = FALSE,
  quote = FALSE, sep = "\t")
