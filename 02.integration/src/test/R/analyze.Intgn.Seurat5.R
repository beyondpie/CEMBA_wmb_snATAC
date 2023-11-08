library(Seurat)
options(Seurat.object.assay.version = "v5")
# Not sure if this is a must
library(SeuratObject)
# Not sure if this is a must
library(SeuratData)
# seurat5 version wrappers
# needed for FastMNNIntegration
library(SeuratWrappers)
# [FIXED] Seurat5 use future inside, this is a must
# for 2 million cells, FindNeighbors needs at least 2G
options(future.globals.maxSize = 5e9)
library(stringr)
library(purrr)
library(ggplot2)
library(Matrix)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  Sa2Integration)
import::from(.from = "utils.R", .directory = rdir,
  setupLogging, closeLogging)

# * configs
projdir <- here::here()
resdir <- file.path(projdir, "19.snap2_integration",
  "out/neuron")
allenS5 <- file.path(resdir, "neuron_10xv3.TF.s5.rds") |> readRDS()
atacS5 <- file.path(resdir, "neuron_atac.TF.s5.rds") |> readRDS()
intgnS5 <- file.path(resdir, "neuron_10xv3_rpca_TF_5.s5.rds") |>
  readRDS()
intgnHS5 <- file.path(resdir, "neuron_10xv3_harmony_TF_5.s5.rds") |>
  readRDS()

# * check if integration finished in rpca
intgnEmbed <- intgnS5@reductions$intgn.rpca@cell.embeddings
allenIntgnEmbed <- intgnEmbed[colnames(allenS5), ]
atacIntgnEmbed <- intgnEmbed[colnames(atacS5), ]
# intgn layer treat first exp as reference, and then map the latter on to
# the same space in RPCAIntegration
# So in this case, allenIntgnEmbed is different with original one
# atacIntgnEmbed is the same as the original one.
sum(allenIntgnEmbed - allenS5@reductions$pca@cell.embeddings)
sum(atacIntgnEmbed - atacS5@reductions$pca@cell.embeddings)

# intgnS5@reductions$pca simply merge the pca results.
sum(intgnS5@reductions$pca@cell.embeddings[colnames(allenS5), ] -
      allenS5@reductions$pca@cell.embeddings)

# * let's follow up for UMAP and clustering
# ** rpca
# now cosuming 21G RAM for 2 million cells with 50 dimensions
# and about one hour
# with default parameter
intgnS5 <- Seurat::RunUMAP(intgnS5,
  reduction = paste0("intgn.", "rpca"), dims = 1:50,
  reduction.name = "umap.rpca")
saveRDS(intgnS5, file.path(resdir, "neuron_10xv3.TF.s5.umap.rds"))

## intgnS5@reductions$umap.rpca = intgnS5@reductions$umap
intgnS5 <- readRDS(file.path(resdir, "neuron_10xv3.TF.s5.umap.rds"))
allenCells <- colnames(intgnS5)[intgnS5@meta.data$modality == "rna"]
atacCells <- colnames(intgnS5)[intgnS5@meta.data$modality == "atac"]

allen.dp <- subset(intgnS5, cells = sample(allenCells, 30000))
atac.dp <- subset(intgnS5, cells = sample(atacCells, 30000))
all.dp <- subset(intgnS5, cells = union(Cells(allen.dp), Cells(atac.dp)))

p.allen <- DimPlot(allen.dp, reduction = "umap.rpca", cols = "red")
p.atac <- DimPlot(atac.dp, reduction = "umap.rpca", cols = "blue")
p.all <- DimPlot(all.dp, reduction = "umap.rpca",
  group.by = "modality")
p <- p.allen + p.atac + p.all
withr::with_pdf("umap_rpca_10xv3_TF.pdf",{
  print(p)
}, width = 30, height = 10)

## saveRDS(object = p, file = "umap_rpca_10xv3_TF_5.rds")
## ggplot2::ggsave(filename = "umap_rpca_10xv3_TF_5.pdf",
  ## plot = p, width = 20, height = 10)
p2 <- DimPlot(intgnS5, reduction = "umap.rpca",
  group.by = "modality",
  combine = TRUE)
ggplot2::ggsave(filename = "umap_rpca_10xv3_TF_5.pdf",
  plot = p2, width = 20, height = 10)

# This step takes at least 30 minutes
# 2 million cells, it takes about 60G RAM during process
intgnS5 <- Seurat::FindNeighbors(intgnS5,
  reduction = paste0("intgn.", "rpca"),  dims = 1:50)

# this step takes about 23G
# it takes at lease 1 hour
intgnS5 <- FindClusters(intgnS5, resolution = 1,
  cluster.name = "rpca_clusters")


# ** harmony
# takes 1.5 hour
intgnHS5 <- Seurat::RunUMAP(
  intgnHS5, reduction = "harmony", dims = 1:50,
  reduction.name = "umap.harmony"
)
saveRDS(intgnHS5, file.path(resdir,
  "neuron_10xv3.Harmony.s5.umap.rds"))

allenCells <- colnames(intgnHS5)[intgnHS5@meta.data$modality == "rna"]
atacCells <- colnames(intgnHS5)[intgnHS5@meta.data$modality == "atac"]

allen.dp <- subset(intgnHS5, cells = sample(allenCells, 30000))
atac.dp <- subset(intgnHS5, cells = sample(atacCells, 30000))
all.dp <- subset(intgnHS5, cells = union(Cells(allen.dp), Cells(atac.dp)))

p.allen <- DimPlot(allen.dp, reduction = "umap.harmony", cols = "red")
p.atac <- DimPlot(atac.dp, reduction = "umap.harmony", cols = "blue")
p.all <- DimPlot(all.dp, reduction = "umap.harmony",
  group.by = "modality")
p <- p.allen + p.atac + p.all
withr::with_pdf("umap_harmony_10xv3_TF.pdf",{
  print(p)
}, width = 30, height = 10)




