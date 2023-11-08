library(tidyverse)
library(ggplot2)
library(ggpubr)
library(fitdistrplus)
packdir <- file.path(here::here(), "package/R")

# * config
projdir <- here::here()
workdir <- file.path(projdir, "04.cCREgene/sa2.cicero")
ciceroDir <- file.path(workdir, "out/sa2pdcsum")

density.theme <-   theme_bw() +
  theme(axis.text.x=element_text(colour = "black", size = 8),
    axis.text.y=element_text(colour = "black", size = 8),
    axis.title.x=element_text(colour = "black", size = 8),
    axis.title.y=element_text(colour = "black", size = 8))


# * functions
loadCor <- function(f) {
  r <- data.table::fread(file = f, header = TRUE, sep = ",",
    data.table = FALSE)
  gene.peak  <- vapply(r$pdc, function(i) {
    t <- strsplit(i, split = "@", fixed = TRUE)
    t[[1]]
  }, c("a", "peak1"))
  r$gene <- gene.peak[1, ]
  r$peak <- gene.peak[2, ]
  return(r)
}

calSingleTailPvalue <- function(x, m, sd)  {
  p_lower <- pnorm(x, mean = m, sd = sd, lower.tail = TRUE, log.p = FALSE)
  p_1_lower <- 1-p_lower
  ind_pos <- (x >= 0)
  p <- ind_pos * p_1_lower + (1-ind_pos) * p_lower
  return(p)
}

# * load data
corRealPear <- loadCor(
  file.path(ciceroDir, "sa2subclass.cor.pearson.real.csv"))
corRdmPear <- loadCor(
  file.path(ciceroDir, "sa2subclass.cor.pearson.rdm.shuf.csv"))


# * plot density
withr::with_pdf(new = file.path(ciceroDir, "sa2.pearson.cor.pdf"),
  code =  {
    data.pear <- rbind(corRealPear, corRdmPear)
    data.pear$class <- c(rep("Real.Pearson", nrow(corRealPear)),
      rep("RdmShuf.Pearson", nrow(corRdmPear)))
    p.pear <- ggplot(data.pear, aes(x = cor, colour = class)) +
      geom_density() +
      ggtitle("Pearson correlation") +
      density.theme
    print(p.pear)
  }
)

# * function
fitRdmCor <- function(corRealFile, corRdmFile, corMethod = "spearman") {
  corReal <- loadCor(corRealFile)
  corReal <- corReal[!is.na(corReal$cor), ]
  corRdm <- loadCor(corRdmFile)
  corRdm <- corRdm[!is.na(corRdm$cor), ]
  fit.rdm.shuf <- fitdistrplus::fitdist(corRdm$cor, "norm", keepdata = FALSE)
  m.rdm.shuf <- fit.rdm.shuf$estimate[[1]]
  sd.rdm.shuf <- fit.rdm.shuf$estimate[[2]]
  para <- data.frame(class = "RdmShufPearson", mean = m.rdm.shuf, sd = sd.rdm.shuf)
  data.table::fwrite(para,
    file.path(ciceroDir,
      paste0("sa2subclass.fitnorm.rdm.shuf.", corMethod,".txt")),
    quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  # * get pvalue and fdr
  corReal$Pval <- calSingleTailPvalue(
    x = corReal$cor, m = m.rdm.shuf, sd = sd.rdm.shuf)
  corReal$FDR <- p.adjust(corReal$Pval, method = "BH")
  data.table::fwrite(
    corReal,
    file = file.path(ciceroDir,
      paste0("sa2subclass.pdc.cor.",corMethod,".fdr.tsv")),
    quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  r.alignv1 <- data.frame(
    conns = gsub("@", "|", corReal$pdc),
    pcc = corReal$cor,
    class = rep("real", nrow(corReal)),
    Pval = corReal$Pval,
    FDR = corReal$FDR
  )

  data.table::fwrite(
    r.alignv1,
    file = file.path(ciceroDir,
      paste0("sa2subclass.pdc.cor.",corMethod, ".fdr.alignv1.tsv")),
    quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}


fitRdmCor(
  corRealFile = file.path(ciceroDir, "sa2subclass.cor.pearson.real.csv"),
  corRdmFile = file.path(ciceroDir, "sa2subclass.cor.pearson.rdm.shuf.csv"),
  corMethod = "pearson")

# * filter pos and neg links
filter.pdc <- function(pdc.with.fdr.file,
                       fdr.thres = 0.01,
                       metaCol = "sa2subclass",
                       corMethod = "pearson",
                       outdir = ciceroDir) {
  r.align1 <- data.table::fread(pdc.with.fdr.file,
    header = TRUE, sep = "\t", quote = "", data.table = FALSE)
  index.pass.fdr <- r.align1$FDR <= fdr.thres
  message(paste(sum(index.pass.fdr), "over", nrow(r.align1), "pdc passes."))
  r <- r.align1[index.pass.fdr, ]
  pos.pdc <- r[r$pcc >= 0, ]
  neg.pdc <- r[r$pcc < 0, ]
  message(paste(nrow(pos.pdc), "pos and", nrow(neg.pdc), "neg."))
  if(!dir.exists(outdir)) {
    dir.create(outdir)
  }
  outFile.pos <- file.path(
    outdir, paste("mba.whole", metaCol, corMethod, "pos.pdc.alignv1.tsv", sep = "."))
  outFile.neg <- file.path(
    outdir, paste("mba.whole", metaCol, corMethod, "neg.pdc.alignv1.tsv", sep = "."))
  data.table::fwrite(pos.pdc, file = outFile.pos, quote = FALSE, sep = "\t",
    col.names = TRUE, row.names = FALSE)
  data.table::fwrite(neg.pdc, file = outFile.neg, quote = FALSE, sep = "\t",
    col.names = TRUE, row.names = FALSE)
  return(list(pos.pdc, neg.pdc))
}

pdc.filter.pearson <- filter.pdc(
  pdc.with.fdr.file = file.path(ciceroDir, "sa2subclass.pdc.cor.pearson.fdr.alignv1.tsv"),
  corMethod = "pearson"
)

# * get pdc that are not pass the filtering.
# this is used for hic background.
r <- data.table::fread(
  file = file.path(ciceroDir, "sa2subclass.pdc.cor.pearson.fdr.alignv1.tsv"),
  header = TRUE, sep = "\t", quote = "", data.table = FALSE)

pdc.rdm <- r[r$FDR > 0.01, ] |>
  arrange(abs(pcc)) 

# # positive pairs
n <- 600000
data.table::fwrite(pdc.rdm[seq_len(n), ],
  file = file.path(ciceroDir, "sa2subclass.pdc.rdm.alignv1.tsv"),
  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")




