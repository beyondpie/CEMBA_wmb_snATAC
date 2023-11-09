library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(future)
library(data.table)

projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env, Sa2PeakCalling)

workdir <- file.path(projdir, "18.snap2_peakcalling")
peaksumdir <- file.path(projdir, "supple.07.peakcalling.allinone")

# * prepare finalized peaks (SPM >= 5)
unionpeaks <- Sa2PeakCalling$readUnionPeakSetFile()
finalpeaks <- subset(unionpeaks, spm >= 5)
fpknms <- rownames(finalpeaks)
fpkbed <- data.frame(
  chrom = finalpeaks$seqnames,
  start = finalpeaks$start,
  end = finalpeaks$end,
  name = fpknms
)
rownames(fpkbed) <- fpkbed$name

write.table(fpkbed, file = file.path(peaksumdir,
  "mba.whole.sa2.final.peak.srt.bed"), quote = FALSE,
  sep = "\t", col.names = FALSE, row.names = FALSE)

# * load finalized peak
fpkbed <- Sa2PeakCalling$loadFinalPeakBed()
fpknms <- fpkbed$name

# * preprare finalized peaks for each pL4 and L4 (binary pmat)
pL4unionpeaks <- Sa2PeakCalling$loadpL4UnionPeakBeds()
groupL4nms <- names(pL4unionpeaks)
cleanpL4nms <- gsub("^.*\\.", "", groupL4nms)
pL4finalpeaks <- lapply(seq_along(pL4unionpeaks), \(i) {
  pL4 <- names(pL4unionpeaks[i])
  rawp <- rownames(pL4unionpeaks[[i]])
  r <- fpknms[fpknms %in% rawp]
  message(length(r), " final peaks in ", pL4)
  return(r)
})
names(pL4finalpeaks) <- groupL4nms

# check we have all the peaks
# all_peaks <- unique(unlist(pL4finalpeaks))

outpL4finalpeak_dir <- file.path(peaksumdir, "sa2.pL4.final.peak.srt")
if (!dir.exists(outpL4finalpeak_dir)) {
  dir.create(outpL4finalpeak_dir)
}

pL4fp_df <- lapply(seq_along(pL4finalpeaks), \(i) {
  rawnm <- names(pL4finalpeaks)[i]
  nm <- cleanpL4nms[i]
  fnm <- file.path(outpL4finalpeak_dir, paste0(nm, ".bed"))
  message(nm, " from ", rawnm, " to outfile: ", fnm)
  rdf <- fpkbed[pL4finalpeaks[[i]], ]
  write.table(rdf, file = fnm, quote = FALSE, sep = "\t",
    col.names = FALSE, row.names = FALSE)
  return(rdf)
})

# add L4 binary pmat
# copy from sa2.Fig2.R in figures
## construct pL4 binary pmat.
pL4s <- unique(atacMeta$pL4)
pL4peakdir <- file.path(projdir,
  "supple.07.peakcalling.allinone", "sa2.pL4.final.peak.srt")
pL4Peaks <- lapply(pL4s, \(pl4) {
  fnm <- file.path(pL4peakdir, paste0(pl4, ".bed"))
  Sa2PeakCalling$readBed4File(fnm)
})
names(pL4Peaks) <- pL4s

bpmat <- vapply(pL4s, \(x) {
  r <- rep(0, nrow(peakBed))
  r[match(pL4Peaks[[x]]$name , peakBed$name)] <- 1
  return(r)
}, rep(0, nrow(peakBed))) |> t()
rownames(bpmat) <- pL4s
colnames(bpmat) <- peakBed$name
saveRDS(bpmat,
  file.path(projdir, "figures", "data", "sa2.binary.pmat.pL4.rds")
)


# * prepare finalized peaks for each subclass (binary pmat)
atacMeta <- readRDS(cembav2env$sa2metaFile)
pL4finalpeaks <- Sa2PeakCalling$loadpL4FinalPeakBeds()
scCol <- "subclass_label_v3"
subclasses <- unique(atacMeta[, scCol])
pL4tosubclass <- unique(atacMeta[, c(scCol, "pL4")])

subclassFinalpeaks <- lapply(seq_along(subclasses), \(i) {
  sc <- subclasses[i]
  pL4s <- pL4tosubclass[pL4tosubclass[, scCol] == sc, 2]
  rawp <- if (length(pL4s) == 1) {
     pL4finalpeaks[[pL4s]]$name
  } else {
    unique(unlist(
      lapply(pL4s, \(x) {
        pL4finalpeaks[[x]]
      })
    ))
  }
  message(sc, " includes: ",
    paste(pL4s, ","), " and has ", length(rawp), " peaks.")
  fpknms[fpknms %in% rawp]
})
names(subclassFinalpeaks) <- subclasses

outscfinalpeak_dir <- file.path(peaksumdir,
  "sa2.subclassv3.final.peak.srt")
if (!dir.exists(outscfinalpeak_dir)){
  dir.create(outscfinalpeak_dir)
}

subclassFinalpeaks_df <- lapply(seq_along(subclasses), \(i) {
  sc <- subclasses[i]
  scnm <- Sa2PeakCalling$renameAllenSubclass(sc)
  fnm <- file.path(outscfinalpeak_dir,
    paste0(scnm, ".bed"))
  message(sc, " to outfile: ", fnm)
  rdf <- fpkbed[subclassFinalpeaks[[sc]], ]
  write.table(rdf, file = fnm, quote = FALSE, sep = "\t",
    col.names = FALSE, row.names = FALSE)
  return(rdf)
})

# * prepare CPM for each pL4 and L4 (pmat with cpm)
## get in python with sa2 and then convert to matrix in R

# * prepare CPM for each subclass (pmat with cpm)
## get in python with sa2 and then convert to matrix in R


# * prepare subclass single-cell level scRNA-seq expression
# downsample to 2000 or 5000 cells per subclass if possible
# use Allen's dataset

# * prepare subclass level scRNA-seq expression
# use Allen's dataset
# under Seuratv5
# get highly variable genes

