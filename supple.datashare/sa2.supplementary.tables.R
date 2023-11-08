source("paper.R")

# * configs
projdir <- here::here()
outdir <- file.path(projdir,
  "supple.datashare",
  "supplementary.tables/October")
fprefix <- "SupplementaryTable"

# * meta
atacMeta <- readRDS(cembav2env$sa2metaFile)

# * Supplementary 3
# See heavySummaryOfTransferLabel.R in for details
# * Supplementary 4
# Annot of Fig1c UMAP from last submission
L1_reso <- cluSumBySa2$resos$L1 
L2_reso <- cluSumBySa2$loadL2Reos()
write.table(L2_reso, file.path(outdir,
  "SupplementaryTable4.L2reso.tsv"),
  sep= "\t",
  col.names = TRUE, row.names = FALSE, quote= FALSE)
L3_reso <- cluSumBySa2$loadL3Resos()
write.table(L3_reso, file.path(outdir,
  "SupplementaryTable4.L3reso.tsv"),
  sep= "\t",
  col.names = TRUE, row.names = FALSE, quote= FALSE)

L4_reso <- cluSumBySa2$loadL4Resos()
write.table(L4_reso, file.path(outdir,
  "SupplementaryTable4.L4reso.tsv"),
  sep= "\t",
  col.names = TRUE, row.names = FALSE, quote= FALSE)

# * Supplementary 7
# under supple.07: sa2.pL4.final.peak.srt
# make it as a binary mat
allpeak <- Sa2PeakCalling$loadFinalPeakBed()
pL4FinalPeaks <- Sa2PeakCalling$loadpL4FinalPeakBeds()
pL4bpmat <- matrix(0, nrow = length(pL4FinalPeaks),
  ncol = nrow(allpeak))
rownames(pL4bpmat) <- names(pL4FinalPeaks)
colnames(pL4bpmat) <- allpeak$name

for(i in seq_along(pL4FinalPeaks)) {
  pl4 <- names(pL4FinalPeaks)[i]
  index <- match(pL4FinalPeaks[[i]]$name, allpeak$name)
  pL4bpmat[pl4, index] <- 1
  message(pl4, ": ", sum(index > 0))
}

# check
sum(colSums(pL4bpmat) > 0)
## # origina one
bpmat.fig2 <- readRDS(
  file.path(projdir, "figures/data", "sa2.binary.pmat.pL4.rds"))
bpmat.fig2 <- bpmat.fig2[rownames(pL4bpmat), colnames(pL4bpmat)]
sum(abs(pL4bpmat-bpmat.fig2))

# 2.9G
write.table(pL4bpmat, file = file.path(
  outdir, "SupplementaryTable7.pL4bpmat.csv"),
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ","
)
# after zip, will be 100M
## library(Matrix)
## bpmat.sparse <- as(pL4bpmat, "sparseMatrix")
## pL4bpmat_outdir <- file.path(outdir, "SupplemenatryTable7")
## if(!dir.exists(pL4bpmat_outdir)) {
##   dir.create(pL4bpmat_outdir)
## }
## writeMM(obj = bpmat.sparse,
##   file = file.path(pL4bpmat_outdir, "pL4_to_peak.sparse.mtx"))
## write.table(rownames(bpmat.sparse),
##   file = file.path(pL4bpmat_outdir, "pL4.rownames.txt"),
##   quote = FALSE, col.names = FALSE, row.names = FALSE)
## write.table(colnames(bpmat.sparse),
##   file = file.path(pL4bpmat_outdir, "peak.colnames.txt"),
##   quote = FALSE, col.names = FALSE, row.names = FALSE)


# * Supplementary 9
pL4tosize <- table(atacMeta$pL4) |>
  toNamedArray.1dtable()
bpmat <- readRDS(
  file.path(projdir,
    "figures", "data", "sa2.binary.pmat.pL4.rds")
)
pL4tonpeak <- rowSums(bpmat)
pL4tosize2npeak <- data.frame(
  pL4 = names(pL4tosize),
  size = pL4tosize,
  npeak = pL4tonpeak[names(pL4tosize)]) |>
  x => x[order(x$size), ]
data.table::fwrite(
  pL4tosize2npeak, file = file.path(outdir, "SuppleTable9.pL4tonpeak.tsv"),
  col.names = TRUE, sep = "\t"
)




# * subclass binary pmat
sctosize <- table(atacMeta$subclass_label_v3) |>
  toNamedArray.1dtable()
subclass2peak <- Sa2PeakCalling$loadSubclassFinalPeakBeds()
sc2npeak <- vapply(subclass2peak, nrow, 1)
names(sc2npeak) <- names(subclass2peak)
names(sc2npeak) <- with(scnmMap,
  subclass_label[match(names(sc2npeak), subclass_label_peak)])
sc2size2npeak <- data.frame(
  subclass_label = names(sctosize),
  size = sctosize,
  npeak = sc2npeak[names(sctosize)]
)  |> x => x[order(x$size), ]

data.table::fwrite(
  sc2size2npeak, file = file.path(outdir, "SuppleTable9.subclass2npeak.tsv"),
  col.names = TRUE, sep = "\t"
)
# then merge into one Excel

# * Supplementary 10
# Module assignments of all the cCREs
statPeakFile <- file.path(
  file.path(projdir, "figures/data", "fig2/sa2_allpeak_nmf",
    "nmfPmat.allpeak.r150.n0.statW"))
statPeak <- loadStatPeak.NMF(statPeakFile)

modOrder <- data.table::fread(
  file = file.path(projdir, "figures/data",
    "fig2", "sa2.allpeak.nmf.module.order.txt"),
  header = FALSE, data.table = FALSE)$V1
mod2rank <- data.frame(
  mod = modOrder,
  rank = seq_len(length(modOrder))
)
rownames(mod2rank) <- mod2rank$mod
peakBed <- Sa2PeakCalling$loadFinalPeakBed()
cCRE2modrank <- data.frame(
  cCRE = statPeak$peak,
  module = mod2rank[statPeak$moduleN ,"rank"]
)

data.table::fwrite(x = cCRE2modrank,
  file = file.path(outdir,
    paste(fprefix, 10, "module.assignments.of.all.cCRE.tsv", sep = "_")),
  sep = "\t",
  col.names = TRUE,
  quote = FALSE)

# * Supplementary 11
# Association of cis regulatory modules with cell subclasses
statSubclass <- loadStatCluster.NMF(
  file = file.path(projdir, "figures/data/fig2",
    "sa2_allpeak_nmf",
    "nmfPmat.allpeak.r150.n0.statH"))
sc2modrank <- data.frame(
  subclass =  statSubclass$cluster,
  module = mod2rank[statSubclass$moduleN, "rank"]
)
data.table::fwrite(x = sc2modrank,
  file = file.path(outdir,
    paste(fprefix, 11, "association.cCREmodule.with.cell.subclass.tsv",
      sep = "_")),
  sep = "\t",
  col.names = TRUE,
  quote = FALSE)


# * Supplementary 12
connStaAll <- data.table::fread(
  file = file.path(projdir, "04.cCREgene/sa2.cicero",
    "out", "sa2pdcsum", "mba.whole.sa2subclass.merge.conns.sta.all"),
  sep = "\t",
  ## NOTE: this is a fixed version and no title
  header = FALSE,
  quote = FALSE,
  data.table = FALSE
) |> setNames(
  object = _,
  nm = c("subclass", "totalconns", "ndistal2distal",
    "nproximal2proximal", "nproximal2distal", "ngeneInpdc", "ncCREsInpdc")
)

data.table::fwrite(
  x = connStaAll,
  file = file.path(outdir, paste0(fprefix, 12, ".tsv")),
  sep = "\t",
  col.names = TRUE
)

# * Supplementary 14
## Module assignment of the putative enhancers
mod.order <- data.table::fread(
  file = file.path(projdir, "figures/data",
    "fig3", "sa2.ppdc.nnf.module.order.txt"),
  header = FALSE,
  data.table = FALSE
)$V1

mod2rank <- data.frame(
  mod = mod.order,
  rank = seq_len(length(mod.order))
)
rownames(mod2rank) <- mod2rank$mod
module <- 54
tag <- "ppdc"
statPeakFile <- file.path(here::here(), "figures/data",
  "fig3/sa2_ppdcpeak_nmf", "nmfPmat.ppdcpeak.r54.n0.statW")
statClusterFile <- file.path(here::here(), "figures/data",
  "fig3/sa2_ppdcpeak_nmf", "nmfPmat.ppdcpeak.r54.n0.statH")
statPeak <- loadStatPeak.NMF(file = statPeakFile)
cCRE2mod <- data.frame(
  cCRE = statPeak$peak,
  mod = mod2rank[statPeak$moduleN , "rank"]
)
data.table::fwrite(
  x = cCRE2mod,
  file = file.path(outdir,
    paste(fprefix, 14,
      "sa2.module.assignment.of.putative.enhancers.tsv", sep = "_")),
  sep = "\t",
  col.names = TRUE
)

# * Supplementary 15
## Module assignment of the subclasses in ppdc
statClusterFile <- file.path(here::here(), "figures/data",
  "fig3/sa2_ppdcpeak_nmf", "nmfPmat.ppdcpeak.r54.n0.statH")
statSubclass <- loadStatCluster.NMF(file = statClusterFile)
sc2mod <- data.frame(
  subclass = statSubclass$cluster,
  mod = mod2rank[statSubclass$moduleN, "rank"]
)
data.table::fwrite(
  x = sc2mod,
  file = file.path(outdir,
    paste(fprefix, 15,
      "module.assignment.of.subclass.in.enhancer.analysis.tsv",
      sep = "_")),
  sep = "\t",
  col.names = TRUE
)

# * Supplementary 16
module <- 54
motifList <- lapply(seq_len(module), function(r) {
  a <- data.table::fread(file = file.path(
   projdir, "figures/data/fig3/nmf.ppdcpeak.r54.motif", 
    paste0("homer_n", r), "knownResults.txt"), sep = "\t",
    header = TRUE, data.table = FALSE)
  motifnms <- vapply(a[,1], function(i) {
    unlist(strsplit(i, split = "/", fixed = TRUE))[1]
  }, FUN.VALUE = "CTCF")
  logp <- a[, 4]
  motifnms <- vapply(a[,1], function(i) {
    unlist(strsplit(i, split = "/", fixed = TRUE))[1]
  }, FUN.VALUE = "CTCF")
  result <- data.frame(motif = motifnms, logp = logp)
  return(result)
})
names(motifList) <- paste0("m", seq_len(module))
all.motif.nms <- unique(unlist(
  lapply(motifList, function(i) {return(i$motif)})
))
cutoff.logp <- (-10)
motifList.filter <- lapply(motifList, function(i) {
  subset(i, logp <= cutoff.logp)
})

out.motif.dir <- file.path(outdir, "outMotif")

invisible(
  lapply(names(motifList.filter), function(i) {
    r <- motifList.filter[[i]]
    j <- mod2rank[i, "rank"]
    data.table::fwrite(x = r,
      file = file.path(out.motif.dir, paste0("nmf.", j, ".tsv")),
      sep = "\t",
      col.names = TRUE
    )
  })
)

# * Supplementary 17
co_outdir <- file.path(outdir, "SupplementaryTable17")
if(!dir.exists(co_outdir)) {
  dir.create(co_outdir)
}
GRNs <- readRDS(
  file.path(projdir, "figures", "data/fig4",
    "sa2.GRNs.baseGRN_subclass.nlogp3_topk10000.rds"))
tags <- map(seq_along(GRNs), \(i) {
  sc <- names(GRNs)[i]
  g <- GRNs[[i]]
  write.table(x = g,
    file = file.path(co_outdir, paste0(sc, ".nlogp3topk10000.csv")),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
  message(sc, " GRN is saved.")
})



# * Supplementary 18
cnt.motif <- readRDS(file.path(projdir, "figures/data/fig4",
  "sa2.all.subclass.cnt.motif.color3.rds"))
write.table(x = cnt.motif, file = file.path(outdir, "SupplementaryTable18.csv"),
  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = ",")
# * Supplementary 19
# meta/sa2.subclass2region2score.csv

# * Supplementary 20
scnmMap <- cembav2env$load.sa2.subclass.names()
rownames(scnmMap) <- scnmMap$subclass_id_label
sc2r <- read.table(
  file = file.path(projdir, "meta", "sa2.subclass2region2score.csv"),
  header = TRUE, sep = ",")
r2color <- S4Vectors::second(largeRegionColors)
names(r2color) <- S4Vectors::first(largeRegionColors)
sc2r$rcolor <- r2color[sc2r$r]
rownames(sc2r)

sc2nt <- unique(atacMeta[,
  c("subclass_label_v3", "NT_v3", "nt_type_color_v3")])
sc2nt$nt_type_color_v3[is.na(sc2nt$nt_type_color_v3)] <- "grey"
rownames(scnmMap) <- scnmMap$subclass_label
sc2nt$subclass <- scnmMap[sc2nt$subclass_label_v3, "subclass_label_peak"]
rownames(sc2nt) <- sc2nt$subclass
ntcolor <- sc2nt[, c("NT_v3", "nt_type_color_v3")] |>
  unique() |>
  x => {
    r <- x[, 2] |>
      setNames(object = _, nm = x[ ,1])
  }

sc2id <- scnmMap[, c("subclass_label_peak", "subclass_id")]
rownames(sc2id) <- sc2id$subclass_label_peak
atacsc_id <- sc2id[rownames(sc2nt), "subclass_id"]
osc <- rownames(sc2nt)[order(atacsc_id)]

cnt <- readRDS(file.path(projdir, "figures/data/fig4",
  "sa2.all.subclass.cnt.motif.color3.rds"))

motri <- df2tri(as.data.frame(t(cnt[, osc[osc %in% colnames(cnt)]])))
colnames(motri) <- c("subclass", "motif", "count")
motri$mainclass <- sc2nt[motri$subclass, "NT_v3"]
motri$region <- sc2r[motri$subclass, "r"]
motri$regionpct <- sc2r[motri$subclass, "r_s"]

data.table::fwrite(
  x = motri,
  file = file.path(
    outdir,
    paste("SupplementaryTable", 19,
      "counts.of.the.network.motifs.in.the mouse.brain.region.tsv", sep = "_")),
  sep = "\t",
  col.names = TRUE
  )

# * 2023-10-31
# * update supplementary table to xlsx as much as possible.
scnmMap <- cembav2env$load.sa2.subclass.names()
scpeak2sc <- scnmMap[ , c("subclass_label_peak", "subclass_label")]
rownames(scpeak2sc) <- scpeak2sc$subclass_label_peak

sc2scpeak <- scpeak2sc
rownames(sc2scpeak) <- sc2scpeak$subclass_label

allenMeta <- Sa2Integration$loadAllenAnnot()
atacMeta <- readRDS(cembav2env$sa2metaFile)
sc2nt <- unique(atacMeta[, c("subclass_label_v3", "NT_v3")])
sc2nt$sc <- sc2scpeak[sc2nt$subclass_label_v3, "subclass_label_peak"]
rownames(sc2nt) <- sc2nt$sc

sc2class <- unique(atacMeta[, c("subclass_label_v3", "class_label_v3")])
sc2class$sc <- sc2scpeak[sc2nt$subclass_label_v3, "subclass_label_peak"]
rownames(sc2class) <- sc2class$sc


library(openxlsx)
# ** Supplementary Table 17
# subclass-specific filtered GRN
GRNs <- readRDS(
  file.path(projdir, "figures", "data/fig4",
    "sa2.GRNs.baseGRN_subclass.nlogp3_topk10000.rds"))
## names(GRNs) <- scpeak2sc[names(GRNs), "subclass_label"]

todf.GRN <- function(GRNs) {
  map(seq_along(GRNs), \(i) {
    sc <- names(GRNs)[i]
    r <- GRNs[[i]]
    ## r$subclass <- scpeak2sc[sc, "subclass_label"]
    r$subclass <- sc
    return(r)
  }) |> x => do.call(rbind, x)
}

# this works in principle, but one nt with > 1 million line
# then open excel will have error.
## NT_v3 <- unique(sc2nt$NT_v3)
## GRNs.byNT <- map(NT_v3, \(nt) {
##   scs <- sc2nt[sc2nt$NT_v3 %in% c(nt), "sc"]
##   todf.GRN(GRNs[intersect(scs, names(GRNs))])
## })
## names(GRNs.byNT) <- NT_v3
## write.xlsx(GRNs.byNT[3],
##   file = file.path(outdir, "Supplementary Table 17.xlsx"))

classes <- unique(sc2class$class_label_v3)
GRNs.byclass <- map(classes, \(nt) {
  scs <- sc2class[sc2class$class_label_v3 %in% c(nt), "sc"]
  todf.GRN(GRNs[intersect(scs, names(GRNs))])
})
names(GRNs.byclass) <- classes
write.xlsx(GRNs.byclass,
  file = file.path(outdir, "Supplementary Table 17.xlsx"))

# * supplementary 16
module <- 54
motifList <- lapply(seq_len(module), function(r) {
  a <- data.table::fread(file = file.path(
   projdir, "figures/data/fig3/nmf.ppdcpeak.r54.motif", 
    paste0("homer_n", r), "knownResults.txt"), sep = "\t",
    header = TRUE, data.table = FALSE)
  motifnms <- vapply(a[,1], function(i) {
    unlist(strsplit(i, split = "/", fixed = TRUE))[1]
  }, FUN.VALUE = "CTCF")
  logp <- a[, 4]
  motifnms <- vapply(a[,1], function(i) {
    unlist(strsplit(i, split = "/", fixed = TRUE))[1]
  }, FUN.VALUE = "CTCF")
  result <- data.frame(motif = motifnms, logp = logp)
  return(result)
})
names(motifList) <- paste0("m", seq_len(module))
all.motif.nms <- unique(unlist(
  lapply(motifList, function(i) {return(i$motif)})
))
## cutoff.logp <- (-10)
## motifList.filter <- lapply(motifList, function(i) {
##   subset(i, logp <= cutoff.logp)
## })

mat.motif <- matrix(data = NA,
  nrow = module, ncol = length(all.motif.nms),
  dimnames = list(paste0("module", seq_len(module)),
    all.motif.nms))

for (i in seq_along(motifList)) {
  mat.motif[i, motifList[[i]]$motif] <- motifList[[i]]$logp
}

mat.motif <- as.data.frame(mat.motif)
write.xlsx(mat.motif,
  file = file.path(outdir, "Supplementary Table 16.xlsx"),
  rowNames = TRUE)
# add README manually (about p-values)

# * supplementary 10
cCRE2modrank <- data.table::fread(
  file = file.path(outdir,
    paste(fprefix, 10, "module.assignments.of.all.cCRE.tsv", sep = "_")),
  sep = "\t",
  header = TRUE,
  data.table = FALSE
)

mouse_chrm <- c(paste0("chr", 1:19), "chrX", "chrY")

cCRE2modrank.bychrom <- map(mouse_chrm, \(chrm) {
  p <- paste0(chrm, ":")
  cCRE2modrank[grepl(p, cCRE2modrank$cCRE),]
})
names(cCRE2modrank.bychrom) <- mouse_chrm
write.xlsx(x = cCRE2modrank.bychrom,
  file = file.path(outdir, "Supplementary Table 10.xlsx"))

# * supplementary 7
cCREs.all <- Sa2PeakCalling$loadFinalPeakBed()

chr_1_10 <- paste0("chr", 1:10, collapse = ":|")
cCREs.chr1_10 <- cCREs.all[grepl(chr_1_10,cCREs.all$name), ]
cCREs.chr10_Y <- cCREs.all[!grepl(chr_1_10,cCREs.all$name), ]
cCREs.list <- list(chr_1_10 = cCREs.chr1_10,
  chr_11_Y = cCREs.chr10_Y)

write.xlsx(x = cCREs.list,
  file = file.path(outdir, "Supplementary Table 7.xlsx"))
