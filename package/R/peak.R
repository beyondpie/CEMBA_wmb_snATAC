#' @export
getPeakCoords <- function(peakdf, sep1 = ":", sep2 = "-") {
  r <- paste(
    peakdf[,1],
    paste(peakdf[,2], peakdf[,3], sep = sep2), sep = sep1)
  return(r)
}

#' @export
normalizePmat <- function(snap,
                          method = "abs",
                          scaleFactor = 1e6) {
  pmat <- snap@pmat
  O <- Matrix::colSums(pmat)
  if(method == "abs") {
    U <- sum(snap@metaData$UQ)
  } else if (method == "rela"){
    U <- sum(O)
  } else {
    stop(method, " is not abs or rela.")
  }
  cpm <- (O / U) * scaleFactor
  return(cpm)
}

#' @export
#' cpm: cell by features
processCPM.log2.scale.cap <- function(cpm,
                       islog2 = TRUE, scaleShape = 1, 
                       isScale = TRUE,
                       isCapped = TRUE,
                       lowCapQuantile = 0.01) {
  postcpm <- cpm
  if(islog2) {
    postcpm <- log2(postcpm + scaleShape)
  }
  if(isScale) {
    # scale for each column, which is opereated on all the rows.
    postcpm <- scale(postcpm, center = TRUE, scale = TRUE)
  }
  if(isCapped) {
    maxVal <- quantile(postcpm, 1-lowCapQuantile)
    minVal <- quantile(postcpm, lowCapQuantile)
    postcpm[postcpm >= maxVal] <- maxVal
    postcpm[postcpm <= minVal] <- minVal
  }
  return(postcpm)
}

#' @export
calCellRatio <- function(snap) {
  binaryPmat <- SnapATAC::makeBinary(obj = snap, mat = "pmat")@pmat
  r <- Matrix::colSums(binaryPmat) / nrow(binaryPmat)
  return(r)
}


#' Return probability that at least half the cells express,
#' if we have observed k of n cells expressing
#' 
#' @param k integer number of cell having a given peak
#' @param n integer number of cell
#' @param f numeric should be [0,1], default is 0.5
#' @param shape1 numeric beta dist shape1, default is 1.5
#' @param shape2 numeric beta dist shape2, default is 2
#'
#' NOTE:
#' p(theta >= 0.5| k, n ), while prior(theta) ~ Beta(shape1, shape2)
#' And we use it to describe the following term:
#' Probability that at least half the cells express,
#' when we observe k positives among n cells.
#' @export
estProbNoLessThanF <- function(k, n,
                              f = 0.5,
                              shape1 = 1.5, shape2 = 2) {
  probLessThanF <- pbeta(f, shape1 + k, shape2 - k + n)
  return(1 - probLessThanF)
}

#' @export
loadParsePmat <- function(celltype,
                          pmatDir,
                          prefix = "mba.whole.celltype.",
                          suffix = ".parsePmat.txt",
                          header = TRUE,
                          cols = c("totalcount", "frac", "cpm", "relanorm")
                          ) {
  f <- file.path(pmatDir, paste0(prefix, celltype, suffix))
  if (!file.exists(f)) {
    warning(celltype, " does not have the pmat yet.")
    return(NULL)
  }
  peakStat <- data.table::fread(
    file = f, header = header, sep = "\t", quote = "",
    data.table = FALSE, verbose = FALSE
  )
  r <- peakStat[, cols, drop = FALSE]
  rownames(r) <- peakStat$peak
  return(r)
}

#' @export
getPmatFromSnap <- function(rawsnap, colnm,
                            grp, nds = NULL, ncore = 1) {
  meta <- snap@metaData
  if(!colnm %in% colnames(snap@metaData)) {
    stop(colnm, " is not in snap.")
  }
  cols <- snap@metaData[, colnm]
  index <- cols %in% grp
  if(!any(index)) {
    stop(grp, " is not in ", colnm)
  }
  
  rows <- if(!is.null(nds)) {
    message("Down sample cells to ", nds)
    sort(sample(x = which(index),
      size = min(sum(index), nds), replace = FALSE))
  } else {
    which(index)
  }
  # * load pmat
  snap <- rawsnap[rows, , drop = FALSE]
  message("Load pmat into snap.")
  snap <- SnapATAC::addPmatToSnap(obj = snap,
    do.par = TRUE, num.cores = ncore)
  message("Load pmat: done.")
  # * get pmat
  r <- snap@pmat
  barcodes <- with(snap@metaData, paste(sample, barcode, sep = "."))
  peaks <- snap@peak$name
  rownames(r) <- barcodes
  colnames(r) <- peaks
  return(r)
}

#' @export
getcov.array <- function(x) {
  cv <- sd(x) / abs(mean(x))
  return(cv)
}


#' Get coefficicent of variation for each features.
#' mat: cell by features
#' @export
getCOV.mat <- function(mat) {
  cov <- apply(mat, 2, getcov.array)
  return(cov)
}

#' @export
loadStatPeak.NMF <- function(file) {
  statPeak <- data.table::fread(
    file = file,
    header = FALSE, sep = "\t", data.table = FALSE
  )
  colnames(statPeak) <- c(
    "peak", "index", "class0", "contributes", "featureScore", "selt_fs_list",
    "selt_med_list")
  statPeak$moduleN <- paste(
    "m", as.integer(statPeak$class0 + 1), sep = "")
  rownames(statPeak) <- statPeak$peak
  return(statPeak)
}

#' @export
loadStatCluster.NMF <- function(file, moduleCol = "class0") {
  statCluster <- data.table::fread(
    file = file, header = TRUE, sep = "\t", data.table = FALSE
  )
  colnames(statCluster) <- c(
    "cluster", "index", "class0", "class1", "contributes",
    "sparseness", "entropy")
  statCluster$moduleN <- paste(
    "m", as.integer(statCluster[, moduleCol] + 1), sep = "")
  rownames(statCluster) <- statCluster$Cluster
  return(statCluster)
}

#' @export
downsample.index.1 <- function(index.all,
                               mod.index,
                               mod.ordered = NULL,
                               score.index = NULL,
                               size = 200,
                               seed = 2022) {
  set.seed(seed = seed)
  mods <- if(!is.null(mod.ordered)) {
    mod.ordered
  } else {
    unique(mod.index)
  }
  index.list <- lapply(mods, function(m) {
    index <- mod.index %in% m
    if(sum(index) < 1) {
      warning("Module ", m, " has no features.")
      return(NULL)
    }
    cols <- which(index)
    cols.sampled <- if(!is.null(score.index)) {
      cols[order(score.index[cols], decreasing = TRUE)][
        1:min(length(cols), size)]
    } else {
      sample(cols,size = min(length(cols), size),
        replace = FALSE)
    }
    return(cols.sampled)
  })
  index.list[sapply(index.list, is.null)] <- NULL
  return(unlist(index.list))
  ## index.sampled <- if( !is.null(mod.ordered) ||
  ##                        !is.null(score.index)) {
  ##   unlist(index.list)
  ## } else {
  ##   sort(unlist(index.list))
  ## }
  ## return(index.sampled)
}

#' @export
downsampleStatPeak.1 <- function(statPeak,
                                 mod.ordered = NULL,
                                 size = 200,
                                 seed = 2022) {
  index.sampled <- downsample.index.1(
    index.all = seq_len(nrow(statPeak)),
    mod.index = statPeak$moduleN,
    mod.ordered = mod.ordered,
    score.index = statPeak$featureScore,
    size = size,
    seed = seed
  )
  statPeak.sampled <- statPeak[index.sampled, ]
  return(statPeak.sampled)
}

#' @export
orderStatPeakByModuleRank <- function(statPeak, mod.ordered){
  index.list <- lapply(mod.ordered, function(m) {
    index.raw <- which(statPeak$moduleN %in% m)
    index <- index.raw[
      order(statPeak$featureScore[index.raw], decreasing = TRUE)]
    return(index)
  })
  index.ordered <- unlist(index.list)
  r <- statPeak[index.ordered, ]
  return(r)
}


#' @export
downsampleMat.column.1 <- function(mat,
                                   mod.col,
                                   score.col= NULL,
                                   mod.ordered = NULL,
                                   size = 200, seed = 2022) {
  index.sampled <- downsample.index.1(
    index.all = seq_len(ncol(mat)),
    mod.index = mod.col,
    mod.ordered = mod.ordered,
    size = size,
    seed = seed,
    score.index = score.col
  )
  mat.sampled <- mat[ , index.sampled]
  mod.col.sampled <- mod.col[index.sampled]
  return(
    list(mat = mat.sampled, modCol = mod.col.sampled)
  )
}

#' @export
orderMatByModuleRank.column <- function(mat,
                                        mod.col,
                                        score.col,
                                        mod.ordered) {
  index.list <- lapply(mod.ordered, function(m) {
    index.raw <- which(mod.col %in% m)
    index <- index.raw[
      order(score.col[index.raw], decreasing = TRUE)]
    return(index)
  })
  index.all <- unlist(index.list)
  return(mat[, index.all])
}

#' @export
getPeakModuleAvgScore.NMF <- function(cpm, moduleOfPeak,
                                      fn = rowMeans) {
  modules <- unique(moduleOfPeak)
  avgScores <- vapply(modules, function(m) {
    a <- cpm[, moduleOfPeak %in% m, drop = FALSE]
    return(fn(a))
  }, FUN.VALUE = rep(0, nrow(cpm)))
  return(avgScores)
}

#' @export
getTopRankClusterForPeakModule <- function(cluster2moduleScore,
                                           modules = colnames(cluster2moduleScore),
                                           topk = 3) {
  n_module <- ncol(cluster2moduleScore)
  r <- vapply(seq_len(n_module), function(i){
    return(
      order(cluster2moduleScore[, i], decreasing = TRUE)[1:topk]
    )
  },FUN.VALUE = rep(1, topk))
  if(topk < 2) {
    r <- matrix(data = r, nrow = 1)
  }
  colnames(r) <- modules
  return(r)
}

#' @export
rankModules <- function(topRank.module,
                        modules = colnames(topRank.module),
                        avg.fn = colMeans) {
  avgRank <- avg.fn(topRank.module)
  index.order <- order(avgRank, decreasing = FALSE)
  return(modules[index.order])
}

#' depends on colFn in utils.R
#' @export
orderModules.cpm.statPeak <- function(cpm, statPeak,
                                      fn.avgModule = rowMeans,
                                      topk = 1,
                                      fn.rank = min) {
  moduleOfPeak <- statPeak$module
  cluster2moduleScore <- getPeakModuleAvgScore.NMF(
    cpm = cpm,
    moduleOfPeak = moduleOfPeak,
    fn = fn.avgModule)
  topRank.module <- getTopRankClusterForPeakModule(
    cluster2moduleScore = cluster2moduleScore,
    modules = colnames(cluster2moduleScore),
    topk = topk
  )
  colMins <- colFn(col.fn = fn.rank)
  modules.order <- rankModules(
    topRank.module,
    modules = colnames(topRank.module),
    avg.fn = colMins)
  return(modules.order)
}

#' Global modules will be ranked at the first positions.
#' @export
orderModules.cpm.statPeak2 <- function(cpm,
                                       statPeak,
                                       fn.avgModule = rowMeans,
                                       topk = 1,
                                       fn.rank = min,
                                       top.quantile = 0.95,
                                       global.cover.ratio = 0.9) {
  moduleOfPeak <- statPeak$module
  cluster2moduleScore <- getPeakModuleAvgScore.NMF(
    cpm = cpm,
    moduleOfPeak = moduleOfPeak,
    fn = fn.avgModule)
  topRank.module <- getTopRankClusterForPeakModule(
    cluster2moduleScore = cluster2moduleScore,
    modules = colnames(cluster2moduleScore),
    topk = topk
  )
  colMins <- colFn(col.fn = fn.rank)
  modules.order <- rankModules(
    topRank.module,
    modules = colnames(topRank.module),
    avg.fn = colMins)
  ## put global modules in front of other modules
  top.value <- quantile(cpm, top.quantile)
  nsubclass.mod <- colSums(cluster2moduleScore >= top.value)
  global.modules <- colnames(cluster2moduleScore)[
    which(nsubclass.mod >= (global.cover.ratio) * nrow(cluster2moduleScore) )
  ]
  mod.reorder <- c(
    global.modules, setdiff(modules.order, global.modules)
  )
  return(mod.reorder)
}

