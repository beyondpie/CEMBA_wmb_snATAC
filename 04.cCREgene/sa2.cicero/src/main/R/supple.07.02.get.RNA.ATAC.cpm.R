library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratObject)

library(Matrix)

projdir <- here::here()
workdir <- file.path(projdir, "04.cCREgene/sa2.cicero")
packdir <- file.path(here::here(), "package", "R")
import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env, Sa2PeakCalling, Sa2Integration)
import::from(.from = "utils.R", .directory = packdir,
  aggMatByCol)
import::from(.from = "annot.R", .directory = packdir,
  mapAllenSubclassTomCGSubclass2)

# * configs
peakfile <- Sa2PeakCalling$finalpeakBedFile
allpdcPairfnm <- cembav2env$sa2.all.pdc.pairs.only

# * load peak info
peakfile <- Sa2PeakCalling$finalpeakBedFile
allPeaks <- data.table::fread(
  file = peakfile, sep = "\t", header = FALSE,
  data.table = FALSE)
colnames(allPeaks) <- c("chr", "start", "end", "name")
peakCoords <- with(allPeaks, paste(chr, paste(start, end, sep = "-"), sep = ":"))
names(peakCoords) <- allPeaks$name

# * load pdc info
all.pdc.pairs <- data.table::fread(
  file = allpdcPairfnm, sep = "\t", header = FALSE,
  data.table = FALSE)
colnames(all.pdc.pairs) <- c("pairs", "proximal")
genes.pdc <- unique(vapply(unique(all.pdc.pairs$pairs),
  function(i) {
    unlist(strsplit(i, split = "|", fixed = TRUE))[1]
  }, FUN.VALUE = "a"))

distal.peaks.pdc <- unique(vapply(unique(all.pdc.pairs$pairs),
  function(i) {
    unlist(strsplit(i, split = "|", fixed = TRUE))[2]
  }, FUN.VALUE = "peak1"))

# * now handle peak cpm data
peakCoords.distal <- peakCoords[distal.peaks.pdc]
cpm.pmat.pbysc <- cembav2env$fastread.csv(
  cembav2env$sa2.cpm.pmat.csv.file)
cpm.pbysc <- cpm.pmat.pbysc[peakCoords.distal, ]

old.rownms <- rownames(cpm.pbysc)
rownames(cpm.pbysc) <- names(peakCoords)[
  match(old.rownms, peakCoords)]
saveRDS(cpm.pbysc,
  file = file.path(workdir, "out/sa2pdcsum",
    "sa2.cpm.pmat.pbysc.distal.rds"))

# * load Allen data
allenseu <- readRDS(cembav2env$sa2.allen.log10cpm.ds1000.seurat)

# * get Allen cpm per subclass
# log(ln, natural log) counts per million
logCPM <- allenseu@assays$RNA$data
avg.logCPM.gbysc <- aggMatByCol(logCPM,
  group.col = allenseu$subclass, fn = rowMeans)
saveRDS(avg.logCPM.gbysc, file = file.path(workdir, "out/sa2pdcsum",
  "sa2.allen.logCPM.gbysc.rds"))

# * select varaibale features from Allen.
allenseu <- Seurat::FindVariableFeatures(allenseu, nfeatures = 5000)
seuvf <- VariableFeatures(allenseu)
allenfeatures <- Sa2Integration$getMarkersList()
length(intersect(seuvf, allenfeatures$k8))
# 283 TF covered
length(intersect(seuvf, allenfeatures$TF))
# 499
length(intersect(allenfeatures$k8, allenfeatures$TF))

allenseu.k8 <- allenseu[allenfeatures$k8, ]
saveRDS(allenseu.k8, file.path(
  here::here(), "22.sa2GRN", "src/main/resource",
  "sa2.allen.logCPM.k8.ds1000.seurat.rds"))

vf.with.tf <- union(
  intersect(seuvf, allenfeatures$k8), allenfeatures$TF)
write.table(vf.with.tf, file.path(here::here(), "22.sa2GRN",
  "src/main/resource", "sa2.allen.vf3281.gene.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

allenseu.vf <- allenseu[vf.with.tf, ]
saveRDS(allenseu.vf, file.path(
  here::here(), "22.sa2GRN", "src/main/resource",
  "sa2.allen.logCPM.vf3281.ds1000.seurat.rds"))

# * save seurat to hdf5
## library(hdf5r)
## conn <- hdf5r::H5File$new(
##   file.path(here::here(), "22.sa2GRN", "src/main/resource",
##     "sa2.allen.logCPM.vf3281.ds1000.cbyg.h5"), mode = "w"
## )
## data.grp <- conn$create_group("X")
#
## data.grp[["logCPM"]] <- allenseu.vf@assays$RNA$data
## conn$close_all()

# * save to mtx like what cell oracle did.
# https://github.com/morris-lab/CellOracle/blob/master/celloracle/data_conversion/seurat_to_mtx.R
# https://github.com/morris-lab/CellOracle/blob/master/celloracle/data_conversion/process_seurat_object.py

library(Matrix)
out_dir <- file.path(here::here(), "22.sa2GRN/src/main/resource",
  "sa2.allen.logCPM.vf3281.ds1000")
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}
# at least 30 minutes, 6.6G
n <- Matrix::writeMM(obj = allenseu.vf@assays$RNA$data,
  file = file.path(out_dir, "gbyc.mtx"))
write.table(colnames(allenseu.vf),
  file = file.path(out_dir, "barcode.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(allenseu.vf),
  file = file.path(out_dir, "gene.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE)

## allenseu.vf <- readRDS(file.path(here::here(),
##   "22.sa2GRN/src/main/resource",
##   "sa2.allen.logCPM.vf3281.ds1000.seurat.rds"))
  
write.table(allenseu.vf@meta.data,
  file = file.path(out_dir, "meta.data.csv"), sep = ",",
  quote = FALSE, row.names = TRUE, col.names = TRUE)

meta <- allenseu.vf@meta.data
meta_data_dtype <- data.frame(row.names = colnames(meta))
for (i in colnames(meta)){
  if (length(class(meta[,i])) == 1){
    class_ = class(meta[,i])
  } else {
    n_class <- length(class(meta[,i]))
    class_ <- class(meta[,i])[n_class]
  }
  meta_data_dtype[i, "dtype"] <- class_
}
write.table(meta_data_dtype,
  file.path(out_dir, "meta_data_dtype.txt"), sep = ",",
  quote = FALSE, row.names = FALSE, col.names = FALSE)

## out_dir <- file.path(here::here(), "22.sa2GRN/src/main/resource",
##   "sa2.allen.logCPM.k8.ds1000")
## if (!dir.exists(out_dir)) {
##   dir.create(out_dir)
## }
## n <- Matrix::writeMM(obj = allenseu.k8@assays$RNA$data,
##   file = file.path(out_dir, "gbyc.mtx"))
## write.csv(colnames(allenseu.k8),
##   file = file.path(out_dir, "barcode.txt"))
## write.csv(rownames(allenseu.k8),
##   file = file.path(out_dir, "gene.txt"))
## write.csv(allenseu.k8@meta.data,
##   file = file.path(out_dir, "meta.data.csv"))

