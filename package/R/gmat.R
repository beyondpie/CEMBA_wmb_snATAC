#' @export
getPseudoCPM1 <- function(gmat,
                    scale.factor = 1e4,
                    asdf = T) {
  totalcnt <- sum(gmat)
  genecnt <- colSums(gmat)
  s <- genecnt  * scale.factor / totalcnt
  names(s) <- colnames(gmat)
  if(asdf) {
    r <- data.frame(c1 = s)
    rownames(r) <- colnames(gmat)
  } else {
    r <- s
  }
  return(r)
}

#' @export
getPseudoCPM <- function(gmat,
                   scale.factor = 1e4,
                   groupby = rep(1, nrow(gmat)),
                   ncores = 1,
                   verbose = FALSE) {
  uGroup <- unique(groupby)
  dfs <- parallel::mclapply(uGroup, function(g) {
    if(verbose) {
      message("Group: ", g)
    }
    mat <- gmat[groupby %in% g, ]
    r <- getPseudoCPM1(gmat = mat, scale.factor = scale.factor, asdf = T)
    return(r)
  }, mc.cores = ncores)
  result <- do.call(cbind, dfs)
  colnames(result) <- as.character(uGroup)
  rownames(result) <- colnames(gmat)
  return(result)
}

#' @export
getPctExp <- function(gmat,
                      groupby,
                      threshold = 0,
                      ncores = 1,
                      verbose = FALSE) {
  uGroup <- unique(groupby)
  dfs <- parallel::mclapply(uGroup, function(g) {
    if(verbose) {
      message("Group: " ,g)
    }
    mat <- gmat[groupby %in% g, ]
    r <- apply(mat, MARGIN = 2, Seurat::PercentAbove, threshold = threshold)
  }, mc.cores = ncores)
  result <- do.call(cbind, dfs)
  colnames(result) <- as.character(uGroup)
  rownames(result) <- colnames(gmat)
  return(result)
}
#' @export
getCPM <- function(gmat,
                   scale.factor = 1e4) {
  r <- gmat * scale.factor / rowSums(gmat)
  return(r)
}

#' @export
loadGmatFromL2 <- function(gmatDir,
                           L2,
                           genes = NULL,
                           dsL3 = -1,
                           dsL2 = -1) {
  files <- list.files(file.path(gmatDir, L2),
    full.names = T, no.. = T)
  if(length(files) < 1) {
    warning("No files found for ", L2)
    return(FALSE)
  }
  snapList <- lapply(files, function(f){
    message("Load file: ", f)
    s <- readRDS(f)
    if(!is.null(genes)) {
      tmp <- s@gmat[, colnames(s@gmat) %in% genes]
      s@gmat <- tmp
    }
    if(dsL3 > 0) {
      index <- sort(sample(
        seq(SnapATAC::nrow(s)), size = min(SnapATAC::nrow(s), dsL3),
        replace = FALSE))
      s <- s[index, , drop = FALSE]
    }
    return(s)
  })
  snap <- SnapATAC::snapListRbind(snapList)
  if(dsL2 > 0){
    index <- sort(sample(
      seq(SnapATAC::nrow(snap)), size = min(SnapATAC::nrow(snap), dsL2),
      replace = FALSE))
    snap <- snap[dsL2, ,drop = FALSE]
  }
  gmat <- snap@gmat
  rownames(gmat) <- paste(snap@sample, snap@barcode, sep = ".")
  return(gmat)
}

#' @export
mergeMatList <- function(matList) {
  nullList <- sapply(matList, is.null)
  if(sum(nullList) == length(matList)) {
    message("All matList is NULL, will return NULL.")
    return(NULL)
  }
  matList <- matList[!nullList]
  mat <- do.call(rbind, matList)
  return(mat)
}

#' @export
getSnapGmat <- function(gmatDir, majorTypePattern = "NonN",
                        barcodes = NULL) {
  # ignore the repeat gmatDir
  majorTypeDirs <- base::dir(path = gmatDir,
    full.names = TRUE, recursive = FALSE, no.. = TRUE)
  matchDirs <- majorTypeDirs[
    grepl(pattern = majorTypePattern, x = majorTypeDirs)]
  matList <- lapply(matchDirs, function(d) {
    message("Current matched dir: ", d)
    subTypeFiles <- list.files(
      path = d, pattern = NULL, full.names = TRUE, no.. = TRUE,
      include.dirs = FALSE
    )
    mList <- lapply(subTypeFiles, function(f) {
      snap <- readRDS(file = f)
      gmat <- snap@gmat
      curBarcodes <- paste(snap@sample, snap@barcode, sep = ".")
      rownames(gmat) <- curBarcodes
      if(!is.null(barcodes)) {
        index <- curBarcodes %in% barcodes
        if(sum(index) > 0 ) {
          gmat <- gmat[index, , drop = FALSE]
        } else {
          gmat <- NULL
        }
      }
      return(gmat)
    })
    mat <- mergeMatList(mList)
    if(is.null(mat)) {
      message("No barcodes selected")
    } else {
      message("Selected barcodes: ", nrow(mat))
    }
    return(mat)
  })
  mat <- mergeMatList(matList)
  if(!is.null(mat) & (!is.null(barcodes))) {
    message("Organize the mat to fit the order of barcodes.")
    curBarcodes <- rownames(mat)
    index <- base::match(x = barcodes, table = curBarcodes, nomatch=0)
    index <- index[index > 0]
    mat <- mat[index, ]
  }
  return(mat)
}

