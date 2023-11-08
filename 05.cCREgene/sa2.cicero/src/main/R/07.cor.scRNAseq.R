library(optparse)
library(Matrix)
pack_dir <- file.path(here::here(), "package/R")
import::from(.from = "utils.R", .directory = pack_dir,
  checkArgsExistOrStop, prepareOutdir, scaleByRow,
  scaleUseRankByRow, fastPearCorOnRow.scaled,
  fastCorOnRow.chunk)
import::from(.from = "annot.R", .directory = pack_dir,
  mapAllenSubclassTomCGSubclass2)
projdir <- here::here()
workdir <- file.path(projdir, "04.cCREgene/sa2.cicero")

op <- list(
  make_option(c("--matRNAFile"), type = "character",
    default = "allen.l2.cpm.ds200.rds"),
  make_option(c("--matATACFile"), type = "character",
    default = "atac.subclass.cpm.rds"),
  make_option(c("--pdcFile"), type = "character",
    default = "pdc.AllenAnnot.rds",
    help = "rdm pdc or real pdc, data.frame with cols of gene,peak"),
  make_option(c("--method"), type = "character",
    default = "pearson"),
  make_option(c("--tag"), type = "character", default = "real"),
  make_option(c("--group"), type = "character",
    default = "AllenAnnotConcat"),
  make_option(c("--shuf"), type = "integer", default = 0),
  make_option(c("--outdir"), type = "character",
    default = "out/AllenAnnotConcat"),
  make_option(c("--chunkSize"), type = "integer", default = 10000)
)

args <- parse_args(OptionParser(option_list = op))

# * manual for pdc
## args$matRNAFile <- file.path(workdir, "out/sa2pdcsum",
##   "sa2.allen.avg.logCPM.gbysc.rds")
## args$matATACFile <- file.path(workdir, "out/sa2pdcsum",
##   "sa2.cpm.pmat.pbysc.distal.rds")
## args$pdcFile <- file.path(workdir, "out/sa2pdcsum",
##   "sa2.g2p.pdc.pair.rds")
## args$method <- "pearson"
## args$tag <- "real"
## args$group <- "subclass"
## args$shuf <- 0
## args$outdir <- file.path(workdir, "out/sa2pdcsum")
## args$chunkSize <- 10000

# * check args
checkArgsExistOrStop(args)

rna.mat.file <- args$matRNAFile
atac.mat.file <- args$matATACFile
pdc.file <- args$pdcFile
outdir <- args$outdir
is.shuf <- ifelse(args$shuf > 0, TRUE, FALSE)
cor.method <- match.arg(
  args$method, choices = c("pearson", "spearman"))

invisible(lapply(c(rna.mat.file, atac.mat.file, pdc.file), function(i){
  if(!file.exists(i)) {
    stop(i, " does not exist.")
  }}))
prepareOutdir(outdir)

tag <- if(is.shuf) {
  paste(args$tag, "shuf", sep = ".")
} else {
  args$tag
}
outfile <- file.path(outdir,
  paste(args$group, "cor", cor.method, tag, "csv",sep = "."))

# * load data
rna.mat <- readRDS(args$matRNAFile)
atac.mat.df <- readRDS(args$matATACFile)
atac.mat <- as.matrix(atac.mat.df)

common.sc <- intersect(colnames(rna.mat), colnames(atac.mat))
# align the column of subclasses
rnaMat <- rna.mat[ , common.sc]
atacMat <- atac.mat[ , common.sc]
# load pdc pairs either in rdm or real
gene2peak <- readRDS(args$pdcFile)

# reduce rnaMat and atacMat
genes <- intersect(
  unique(gene2peak[,1]), rownames(rnaMat))
if(!all(gene2peak[,2] %in% rownames(atacMat))) {
  stop("Some peaks doest not have infomation on atacMat..")
}

rnaMat <- rnaMat[genes, ]
gene2peak <- gene2peak[(gene2peak[,1] %in% rownames(rnaMat)) , ]

# * test
## row.index <- sample(
##   min(nrow(rnaMat),nrow(atacMat)),
##   args$chunkSize)
## p12pairs <- gene2peak[row.index, ]
## cor.method <- "spearman"

# * normalize mat
message("log transformation")
message("current Allen is logCPM, so keep it.")
rnaMat <- rnaMat
atacMat <- log1p(atacMat)
if (is.shuf) {
  message("Will shuf the columns of atacMat")
  set.seed(2023)
  origCols <- colnames(atacMat)
  atacMat <- atacMat[, sample(ncol(atacMat))]
  colnames(atacMat) <- origCols
}

message("Normalizing rna mat.")
s.rna.mat <- if(cor.method == "pearson") {
  scaleByRow(mat = rnaMat)
} else {
  scaleUseRankByRow(mat = rnaMat)
}

message("Normalizing atac mat.")
s.atac.mat <- if(cor.method == "pearson") {
  scaleByRow(mat = atacMat)
} else {
  scaleUseRankByRow(mat = atacMat)
}

# * get correlation
r <- fastCorOnRow.chunk(
  p1byc = s.rna.mat, p2byc = s.atac.mat, p12pairs = gene2peak,
  chunkSize = args$chunkSize, corMethod = cor.method)

# * test real correlation
## r_real <- vapply(seq_len(nrow(p12pairs)), function(i) {
##   gene <- p12pairs$gene[i]
##   peak <- p12pairs$peak[i]
##   cor(s.rna.mat[gene, ], s.atac.mat[peak,], method = cor.method)
## }, FUN.VALUE = 0.1)
## names(r_real) <- with(p12pairs, paste(gene, peak, sep = "@"))
## plot(r, r_real)

# * save result
rdf <- data.frame(
  pdc = names(r),
  cor = r
)
write.table(x = rdf, file = outfile,
  sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Done")
