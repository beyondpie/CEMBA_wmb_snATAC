# This script is prepared by Sai.
library(Matrix)
library(monocle3)
# should be cicero for monocle3
# devtools::install_github(
#  "cole-trapnell-lab/cicero-release", ref = "monocle3")
library(cicero)

projdir <- here::here()
rscdir <- file.path(projdir, "src/main/resource")

args <- commandArgs(trailingOnly = TRUE)
CellType <- args[1]
print(CellType)
path <- "/oasis/tscc/scratch/smamde/mouse_atlas/mtx_files/"
# Read in matrix data using the Matrix package
indata <- Matrix::readMM(paste0(path,CellType,".mtx"))
# Binarize the matrix
indata@x[indata@x > 0] <- 1
indata <- t(indata)
# Format cell info
cellinfo <- read.table(paste0(path,CellType,"_Barcodes.tsv"))
row.names(cellinfo) <- cellinfo$X0
names(cellinfo) <- "cells"

# Format peak info
peakinfo <- read.table("/oasis/tscc/scratch/smamde/mouse_atlas/peaks.tsv")
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
cell_metadata = cellinfo,
gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP


cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

chromosome_length <- read.table("/oasis/tscc/scratch/smamde/mouse_atlas/mm10_chromosome_length.txt")

conns <- run_cicero(cicero_cds, chromosome_length)


saveRDS(conns,paste0(path,CellType,"_cicero_connections.Rds"))

all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = paste0(path,CellType,"_all_peaks.csv"))
write.csv(x = conns, file = paste0(path,CellType,"_cicero_connections.csv"))

