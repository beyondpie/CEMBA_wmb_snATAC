#' Normalization single-nuclei methylation gene-score data (gmat).
#'
#' We will use a transform function on the gmat, then
#' log10CPM normalizaze (option, default not run), scale,
#' findvariablefeatures, and get pca with all the features.
#'
#' @return Seurat object
#' @import Seurat
#' @export
normalizeMethyGmat <- function(snmcSeurat,
                               runCPM = FALSE,
                               rawAssay = "snmcGmat",
                               rawSlot = "counts",
                               doTransform = TRUE,
                               transform = function(x) {
                                 -x
                               },
                               newAssay = "negGmat",
                               nFeature = 2000,
                               vfMethod = "disp",
                               ...) {
  if (doTransform) {
    rawGmat <- Seurat::GetAssayData(
      object = snmcSeurat,
      slot = rawSlot,
      assay = rawAssay
    )
    newGmat <- transform(rawGmat)
    seuratObj <- snmcSeurat
    seuratObj[[newAssay]] <- Seurat::CreateAssayObject(
      counts = newGmat
    )
    Seurat::DefaultAssay(seuratObj) <- newAssay
  } else {
    seuratObj <- snmcSeurat
  }
  if (runCPM) {
    seuratObj <- Seurat::NormalizeData(
      object = seuratObj,
      ...
    )
  }
  seuratObj <- Seurat::ScaleData(
    object = seuratObj,
    rownames(seuratObj),
    ...
  )
  seuratObj <- Seurat::FindVariableFeatures(
    object = seuratObj,
    selection.method = vfMethod,
    nfeatures = nFeature,
    ...
  )
  seuratObj <- Seurat::RunPCA(
    seuratObj,
    features = Seurat::VariableFeatures(seuratObj),
    ...
  )
  return(seuratObj)
}

#' @return data.table
#' @import data.table
#' @export
downsampleSnmcGmat <- function(gmat, meta, l = "L2", n = 200) {
  if(data.table::key(gmat) != "cell") {
    data.table::setkey(gmat, cell)
  }
  a <- meta[, .(cell = sample(cell,
    size = min(length(cell), n),
    replace = F
  )), by = l]
  gmat <- gmat[a$cell, ]
  ## gmat[, cell:=NULL]
  return(gmat)
}

#' @return data.table, cell by gene
mapGmatToGeneName <- function(gmat,
                              gencode,
                              cells,
                              fn = base::rowMeans) {
  geneIds <- colnames(gmat)
  geneNames <- gencode[geneIds, gene_name]
  dupIndex <- duplicated(geneNames)
  dupNames <- unique(geneNames[dupIndex])

  gdupl <- lapply(dupNames, function(g) {
    v <- gmat[, .(fn(do.call(cbind, .SD))),
      .SDcols = which(geneNames %in% g)
    ]
    colnames(v) <- g
    return(v)
  })
  gdup <- do.call(cbind, gdupl)

  gdup$cell <- cells
  data.table::setkey(gdup, cell)
  guniq <- gmat[, .SD, .SDcols = which(!(geneNames %in% dupNames))]
  colnames(guniq) <- geneNames[!(geneNames %in% dupNames)]
  guniq$cell <- cells
  data.table::setkey(guniq, cell)
  r <- guniq[gdup, on = .(cell = cell)]
  r[, cell := NULL]
  return(r)
}

#' @param gmat data.table
#' @param meta data.table
#' @param gencode data.table
#' @return Seurat object with down-sampling
#' @import data.table
#' @export
getMethySeurat <- function(gmat,
                           meta,
                           gencode,
                           n = 200,
                           l = "L2",
                           assay = "snmcGmat",
                           clusterColName = "cluster") {
  message("Down sampling...")
  gmatds <- downsampleSnmcGmat(
    gmat = gmat,
    meta = meta, l = l, n = n
  )
  message("Convert to GeneId to GeneName...")
  cells <- gmatds$cell
  gmatds[, cell := NULL]
  gmatdss <- mapGmatToGeneName(
    gmat = gmatds,
    gencode = gencode,
    cells = cells,
    fn = base::rowMeans
  )

  message("To Seurat ...")
  clusters <- as.data.frame(meta[cells, l, with = FALSE])
  t <- t(gmatdss)
  colnames(t) <- cells
  seurat <- Seurat::CreateSeuratObject(
    counts = t,
    assay = assay
  )
  seurat <- Seurat::AddMetaData(
    object = seurat,
    metadata = clusters[, 1],
    col.name = clusterColName
  )
  return(seurat)
}

#' Integrate single-nuclei ATAC-seq data with single-nuclei
#' methylation-seq data.
#' @param anchorNormMethod characters, "LogNormalize" (default) or "SCT"
#' @return Seurat object
#' @import Seurat
#' @export
integrateSnatacWithSnmc <- function(seuratList,
                                    expSlot = "experiment",
                                    expName = c("mc", "atac"),
                                    defaultAssay = c(NULL, NULL),
                                    preprocess = c(T, T),
                                    anchorNormMethod = "LogNormalize",
                                    pcaDims = 1:50,
                                    anchorDims = 1:50,
                                    kAnchor = 5,
                                    kFilter = 200,
                                    kScore = 30,
                                    ref = NULL,
                                    atacUsePCA = FALSE,
                                    reso = 0.5,
                                    anchorReduction = "cca",
                                    useFeature = "atac",
                                    orFeatures = NULL,
                                    nFeature = 2000,
                                    integrateAssayName = "coSeurat") {
  for (i in 1:length(seuratList)) {
    if (!is.null(defaultAssay[i])) {
      Seurat::DefaultAssay(seuratList[[i]]) <- defaultAssay[i]
    }
    seuratList[[i]][[expSlot]] <- expName[i]
  }

  # * Preprocess data
  for (i in 1:length(seuratList)) {
    if (!preprocess[i]) {
      next
    }
    message(
      "Perfoming data preprocessing for: ",
      names(seuratList[i])
    )
    seuratList[[i]] <- Seurat::NormalizeData(seuratList[[i]])
    seuratList[[i]] <- Seurat::ScaleData(
      seuratList[[i]],
      rownames(seuratList[[i]])
    )
    seuratList[[i]] <- Seurat::FindVariableFeatures(
      seuratList[[i]],
      selection.method = "vst",
      nfeatures = nFeature
    )
    if (names(seuratList[i]) == "atac" & (!atacUsePCA)) {
      next
    } else {
      seuratList[[i]] <- Seurat::RunPCA(
        seuratList[[i]],
        features = rownames(seuratList[[i]])
      )
    }
  }
  # * integration
  if (useFeature %in% names(seuratList)) {
    message("Use features from: ", useFeature)
    features <- Seurat::VariableFeatures(seuratList[[useFeature]])
  } else if (useFeature %in% "rna") {
    message("Use features from RNA")
    features <- orFeatures
  } else {
    message("Use select integration features.")
    features <- Seurat::SelectIntegrationFeatures(seuratList,
      nfeatures = nFeature
    )
  }
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = seuratList,
    dims = anchorDims,
    anchor.features = features,
    normalization.method = anchorNormMethod,
    reduction = anchorReduction,
    reference = ref,
    scale = FALSE,
    k.score = kScore,
    k.filter = kFilter,
    k.anchor = kAnchor
  )
  coSeurat <- Seurat::IntegrateData(
    anchorset = anchors,
    new.assay.name =
      integrateAssayName,
    dims = pcaDims,
    normalization.method =
      anchorNormMethod,
    features = features,
    features.to.integrate = features
  )
  Seurat::DefaultAssay(coSeurat) <- integrateAssayName
  if (anchorNormMethod != "SCT") {
    coSeurat <- Seurat::ScaleData(
      coSeurat,
      features
    )
  }
  coSeurat <- Seurat::RunPCA(coSeurat,
    npcs = max(pcaDims),
    features = features
  )
  coSeurat <- Seurat::RunUMAP(coSeurat,
    reduction = "pca",
    dims = pcaDims
  )
  coSeurat <- Seurat::FindNeighbors(coSeurat, dims = pcaDims)
  coSeurat <- Seurat::FindClusters(coSeurat, resolution = reso)
  return(coSeurat)
}


integrateWithScRNASeq <- function(snapSeurat,
                                  rnaSeurat,
                                  eigDims = NULL,
                                  snapAssay = "GeneScore",
                                  preprocessSnap = TRUE,
                                  preprocessRNA = TRUE,
                                  rnaTypeColnm = "ClusterName",
                                  reso = 0.5,
                                  reduction = "cca",
                                  rnaRunPCA = T,
                                  atacRunPCA = F) {
  if (nrow(snapSeurat) < 1) {
    stop("No cells in snapSeurat.")
  }
  if (nrow(rnaSeurat) < 1) {
    stop("No cells in rnaSeurat.")
  }
  Seurat::DefaultAssay(snapSeurat) <- snapAssay
  if (is.null(eigDims)) {
    eigDims <- 1:nrow(snapSeurat[["pca"]])
  }
  if (preprocessSnap) {
    snapSeurat <- NormalizeData(snapSeurat)
    snapSeurat <- FindVariableFeatures(snapSeurat)
    snapSeurat <- ScaleData(snapSeurat,
      features = rownames(snapSeurat)
    )
    if (atacRunPCA) {
      snapSeurat <- RunPCA(snapSeurat,
        features = VariableFeatures(snapSeurat)
      )
    }
  }
  if (preprocessRNA) {
    rnaSeurat <- NormalizeData(rnaSeurat)
    rnaSeurat <- FindVariableFeatures(rnaSeurat)
    rnaSeurat <- ScaleData(rnaSeurat, features = rownames(rnaSeurat))
    if (rnaRunPCA) {
      rnaSeurat <- RunPCA(rnaSeurat,
        features =
          VariableFeatures(rnaSeurat)
      )
    }
  }
  anchors <- FindTransferAnchors(
    reference = rnaSeurat,
    query = snapSeurat,
    features = VariableFeatures(rnaSeurat),
    reference.assay = "RNA",
    query.assay = snapAssay,
    reduction = reduction
  )
  ## predict type
  transferLabel <- TransferData(
    anchorset = anchors,
    refdata = rnaSeurat@meta.data[, rnaTypeColnm],
    weight.reduction = snapSeurat[["pca"]],
    dims = eigDims
  )
  rownames(transferLabel) <- colnames(snapSeurat)
  t <- data.frame(
    predictId = transferLabel$predicted.id,
    predictMaxScore = apply(transferLabel[, -1], 1, max),
    row.names = colnames(snapSeurat)
  )
  snapSeurat <- AddMetaData(snapSeurat, metadata = t)
  ## impute gene expression
  refdata <- GetAssayData(
    object = rnaSeurat,
    assay = "RNA",
    slot = "data"
  )
  geneImpute <- TransferData(
    anchorset = anchors,
    refdata = refdata,
    weight.reduction = snapSeurat[["pca"]],
    dims = eigDims
  )
  snapSeurat[["RNA"]] <- CreateAssayObject(counts = geneImpute@data)
  rm(geneImpute)

  ## coEmbed
  coEmbed <- merge(x = snapSeurat, y = rnaSeurat)
  DefaultAssay(coEmbed) <- "RNA"
  coEmbed$tech <- ifelse(!is.na(coEmbed@meta.data[, rnaTypeColnm]), "RNA", "ATAC")
  coEmbed <- ScaleData(coEmbed,
    features = VariableFeatures(rnaSeurat),
    do.scale = FALSE
  )
  coEmbed <- RunPCA(coEmbed,
    features = VariableFeatures(rnaSeurat),
    verbose = FALSE
  )
  coEmbed <- RunUMAP(coEmbed, dims = eigDims)
  coEmbed <- FindNeighbors(coEmbed, dims = eigDims)
  coEmbed <- FindClusters(coEmbed, resolution = reso)
  return(coEmbed)
}
setHeuristicColnms <- function(m) {
  rownmMaxOlv <- sapply(colnames(m), function(i) {
    names(which.max(m[, i]))
  })
  colnm <- unlist(lapply(rownames(m), function(i) {
    names(rownmMaxOlv[rownmMaxOlv == i])
  }))
  return(colnm)
}

downsampleSnap <- function(s, ndp = 200) {
  n <- nrow(s)
  if(n <= ndp) {
    return(s)
  }
  return(s[sample(seq(n), size = ndp, replace = FALSE), , drop = FALSE])
}

downsampleSnapList <- function(snapList, ndp = 200) {
  sl <- lapply(snapList, function(s) {
    n <- nrow(s)
    if(n <= ndp) {
      return(s)
    }
    return(s[sample(seq(n), size = ndp, replace = FALSE), , drop = FALSE])
  })
  sR <- SnapATAC::snapListRbind(sl)
  sR <- pureRUtils::orderSnap(sR)
  return(sR)
}

loadNamedMatFromCSV <- function(matFile) {
  t <- data.table::fread(matFile)
  geneNames <- colnames(t)[-1]
  cellNames <- t$name
  mat <- as.matrix(t[,-1])
  rownames(mat) <- cellNames
  colnames(mat) <- geneNames
  return(mat)
}

getGeneSymbol <- function(snapGeneSymbols, rnaGeneEnsmus, ensmu2geneSymbol){
  rnaGeneSymbols <- ensmu2geneSymbol[rnaGeneEnsmus]
  return(intersect(snapGeneSymbols, rnaGeneSymbols))
}

#' @export
get.downsample.fun <- function(minNum = 500,
                               maxNum = 1000){
  f <- function(index, labels) {
    stat <- table(labels)
    ## put divided in the first since sometimes the full production may
    ## introduce overflow problem in R since
    # only 32-bit integer are supported
    idealNum <- (stat / sum(stat)) * minNum * length(stat)
    nms <- names(idealNum)
    lindex <- lapply(nms, function(i) {
      curIndex <- index[which(labels %in% i)]
      if ( stat[i] <= minNum) {
        return(curIndex)
      } else {
        tmpIndex <- sample(curIndex, size = idealNum[i], replace = FALSE)
        tmp2Index <- tmpIndex
        if(length(tmpIndex) <= minNum) {
          tmp2Index <- sample(curIndex, size = minNum, replace = FALSE)
        }
        if (length(tmpIndex) > maxNum) {
          tmp2Index <- sample(tmpIndex, size = maxNum, replace = FALSE)
        }
        return(tmp2Index)
      }
    })
    return(sort(unlist(lindex)))
  }
  return(f)
}

#' Convert snap object to Seurat by gmat
#'
#' @param snap snap object defined in SnapATAC package
#' snap@gmat must be not empty; snap@smat@dmat should not be empty
#' @param eigDims vector, used for choosing PCA components, default 1:50
#' @param assay characters, name used in Seurat object
#' @param pcaPrefix characters, default "SnapATAC_"
#' @param nameDelim characters, default "-"
#' @param useSnapATACEmbed bool, use SnapATAC embedding as pca
#' embeding in Seurat, default TRUE.
#' @return Seurat object
#' @import Seurat
#' @export
snapGmat2Seurat <- function(snap, eigDims = 1:50,
                            assay = "GeneScore",
                            pcaPrefix = "SnapATAC_",
                            nameDelim = "-",
                            useSnapATACEmbed = TRUE) {
  # check snap@gmat
  # check snap@smat@dmat
  metaData <- snap@metaData
  rownames(metaData) <- paste(metaData$sample, metaData$barcode, sep = ".")
  gmatUse <- Matrix::t(snap@gmat)
  colnames(gmatUse) <- paste(metaData$sample, metaData$barcode, sep = ".")
  snapSeurat <- Seurat::CreateSeuratObject(counts = gmatUse, assay =
                                                               assay,
                                           names.delim = nameDelim)
  snapSeurat <- Seurat::AddMetaData(snapSeurat, metadata = metaData)
  if (useSnapATACEmbed) {
    message("Use SnapATAC Embed as pca in Seurat.")
    pcaUse <- snap@smat@dmat[, eigDims]
    rownames(pcaUse) <- paste(metaData$sample, metaData$barcode, sep = ".")
    colnames(pcaUse) <- paste0(pcaPrefix, 1:ncol(pcaUse))
    snapSeurat[["pca"]] <- methods::new(
      Class = "DimReduc", cell.embeddings = pcaUse,
      feature.loadings = matrix(0, 0, 0),
      feature.loadings.projected = matrix(0, 0, 0),
      assay.used = assay, stdev = rep(1, ncol(pcaUse)),
      key = pcaPrefix,
      jackstraw = new(Class = "JackStrawData"), misc = list())
  }
  return(snapSeurat)
}

#' Calculate Overlap Score matrix for ATAC and RNA jointly embedding metaTable
#' @description
#' Based on the co-embeding clustering result, we summation the minmum scores of
#' the percentages of cells on either ATAC or RNA group in each cluster.
#' @param meta data.frame,
#' three columns as ident, atacCluster, rnaCluster defined as in the following
#' parameters
#' @param ident characters, name of column ident, i.e. the cluster Ids for
#' the co-embedded seurat, "coembed.idents" as default.
#' @param atacCol characters, name of column atac,
#' "MajorType"as default
#' @param rnaCol characters, name of column rna,
#' "ClusterName" as default
#' @return data.frame of numeric, atac by rna
#' @export
getOverlapMatrix <- function(meta,
                             ident = "coembed.idents",
                             atacCol = "MajorType",
                             rnaCol = "ClusterName") {
  ident2rna <- data.frame(idents = meta[[ident]], rna_label = meta[[rnaCol]])
  ident2rna <- ident2rna[stats::complete.cases(ident2rna), ]
  ident2atac <- data.frame(idents = meta[[ident]], atac_label = meta[[atacCol]])
  ident2atac <- ident2atac[stats::complete.cases(ident2atac), ]
  rnaTable <- table(ident2rna)
  atacTable <- table(ident2atac)
  rnaPct <- apply(rnaTable, 2, function(x) {
    x / sum(x)
  })
  atacPct <- apply(atacTable, 2, function(x) {
    x / sum(x)
  })
  rnaClusterName <- colnames(rnaPct)
  atacClusterName <- colnames(atacPct)
  calOvlpScoreElement <- function(t1l, t2l) {
    t1PctDF <- data.frame(rnaPct[, t1l])
    colnames(t1PctDF) <- "t1"
    t1PctDF$ident <- rownames(t1PctDF)
    t2PctDF <- data.frame(atacPct[, t2l])
    colnames(t2PctDF) <- "t2"
    t2PctDF$ident <- rownames(t2PctDF)
    comp <- plyr::join(t1PctDF, t2PctDF, by = "ident", type = "full")
    comp[is.na(comp)] <- 0
    comp$ident <- NULL
    comp <- t(comp)
    return(sum(apply(comp, 2, min)))
  }
  rna2atacType <- outer(rnaClusterName, atacClusterName,
                        FUN = paste, sep = "|")
  ovlpScore <- apply(rna2atacType, MARGIN = c(1, 2), FUN = function(i) {
    t <- unlist(strsplit(i, split = "|", fixed = TRUE))
    calOvlpScoreElement(t[1], t[2])
  })
  rownames(ovlpScore) <- rnaClusterName
  colnames(ovlpScore) <- atacClusterName
  return(t(ovlpScore))
}

#' @return data.table, cell by gene
#' @export
mapGmatToGeneName <- function(gmat,
                              gencode,
                              cells,
                              fn = base::rowMeans) {
  geneIds <- colnames(gmat)
  geneNames <- gencode[geneIds, gene_name]
  dupIndex <- duplicated(geneNames)
  dupNames <- unique(geneNames[dupIndex])

  gdupl <- lapply(dupNames, function(g) {
    v <- gmat[, .(fn(do.call(cbind, .SD))),
      .SDcols = which(geneNames %in% g)
    ]
    colnames(v) <- g
    return(v)
  })
  gdup <- do.call(cbind, gdupl)

  gdup$cell <- cells
  data.table::setkey(gdup, cell)
  guniq <- gmat[, .SD, .SDcols = which(!(geneNames %in% dupNames))]
  colnames(guniq) <- geneNames[!(geneNames %in% dupNames)]
  guniq$cell <- cells
  data.table::setkey(guniq, cell)
  r <- guniq[gdup, on = .(cell = cell)]
  r[, cell := NULL]
  return(r)
}

#' @return data.table
#' @import data.table
#' @export
downsampleSnmcGmat <- function(gmat, meta, l = "L2", n = 200) {
  if(data.table::key(gmat) != "cell") {
    data.table::setkey(gmat, cell)
  }
  a <- meta[, .(cell = sample(cell,
    size = min(length(cell), n),
    replace = F
  )), by = l]
  gmat <- gmat[a$cell, ]
  ## gmat[, cell:=NULL]
  return(gmat)
}

#' @export
downSampleOnSnap <- function(snap, cluster, n = 200) {
  snapList <- lapply(unique(cluster), function(i) {
    s <- snap[cluster %in% i, , drop = FALSE]
    if (nrow(s) <= n) {
      return(s)
    } else {
      return(s[sample(nrow(s), size = n, replace = F), ])
    }
  })
  return(SnapATAC::snapListRbind(snapList = snapList))
}

#' @export
getAllenSubclassGeneMarkers <- function(allenMeta,
                                        all.genes,
                                        subclass = "LA-BLA-BMA-PA Glut") {
  index.allen <- allenMeta$upl2 %in% subclass
  # gene from name
  name.pool <- as.vector(
    stringr::str_split(subclass, pattern = " ", simplify = TRUE))
  name.index <- name.pool %in% all.genes
  name.markers <- if (sum(name.index) > 0) {
    message(subclass, " has match genes in its name.")
    name.pool[name.index]
  } else {
    message(subclass, " has no matched gene in its name.")
    NULL
  }

  # top.markers
  top.markers.pool <- as.vector(stringr::str_split(
    allenMeta$top.markers[index.allen], pattern = ",",
    simplify = TRUE))
  top.markers.pool <- top.markers.pool[nchar(top.markers.pool) > 0]
  top.markers.pool <- top.markers.pool[top.markers.pool %in% all.genes]
  top.markers <- names(sort(table(top.markers.pool), decreasing = TRUE))
  message(subclass, " has ", length(top.markers), " top markers.")

  # combo markers
  combo.markers.pool <- as.vector(stringr::str_split(
    allenMeta$combo.markers[index.allen], pattern = ",",
    simplify = TRUE))
  combo.markers.pool <- combo.markers.pool[
    nchar(combo.markers.pool) > 0]
  combo.markers.pool <- combo.markers.pool[
    combo.markers.pool %in% all.genes]
  combo.markers <- names(sort(table(combo.markers.pool),
    decreasing = TRUE))
  message(subclass, " has ", length(combo.markers), " combo markers.")
  # combine all the signals
  signal.nms <- unique(c(name.markers, top.markers, combo.markers))
  return(signal.nms)
}

#' @export
getIntOvlpMat <- function(class = "nn",
                          k.anchor = 5,
                          k.feature = 200,
                          intEmbed = "cca",
                          intRef = "notatac",
                          featuresFrom = "notatac",
                          atacEmbed = "pca",
                          annot1 = "L3",
                          annot2 = "l2",
                          r = 0.8) {
  rdir <- file.path(here::here(), "03.integrate",
    paste("int", class, "L3", sep = "."), "out")
  fnm <- stringr::str_glue(
    "kAnchor{k.anchor}-kFeature{k.feature}",
    "{intEmbed}_Int-ref_{intRef}",
    "feat_{featuresFrom}",
    "atacEmbed_{atacEmbed}",
    "annot-{annot1}_{annot2}.rds",
    .sep = "-")
  # load data
  result <- readRDS(file.path(rdir, fnm))
  mat <- result$intgnClusterList[[paste0("r", r)]][["ovlpMat"]]
  return(mat)
}

