#' @export
splitByChunk <- function(nondecVec, chunkSize = 2012) {
  n <- ceiling(max(nondecVec) / chunkSize)
  r <- lapply(seq(n), function(i) {
    left <- chunkSize * (i-1)
    right <- chunkSize *i
    index <- (nondecVec > left) & (nondecVec <= right)
    if(sum(index) < 1) {
      return(NULL)
    } else {
      return(nondecVec[index])
    }
  })
  rr <- r[!sapply(r, is.null)]
  return(rr)
}

#' @export
h5adToMat.Sparse <- function(filenm) {
  conn <- hdf5r::H5File$new(filename = filenm,
    mode = "r")
  data <- conn[["X"]][["data"]][]
  indices <- conn[["X"]][["indices"]][]
  indptr <- conn[["X"]][["indptr"]][]
  obs <- conn[["obs"]][["_index"]][]
  var <- conn[["var"]][["_index"]][]
  mat <- Matrix::sparseMatrix(i = indices,
    p = indptr, x = data, index1 = FALSE,
    dims = c(length(obs), length(var)))
  rownames(mat) <- obs
  colnames(mat) <- var
  return(mat)
}

#' @export
subset.snap.gmat.h5ad <- function(barcodes, h5adDir, ncore = 2,
                                  is.parallel = FALSE) {
  sampleFiles <- list.files(path = h5adDir,
    full.names = TRUE, no.. = TRUE)
  tmp <- h5adToMat.Sparse(filenm = sampleFiles[1])
  genes <- colnames(tmp)
  rm(tmp)
  
  matList <- if (is.parallel) {
    parallel::mclapply(
      sampleFiles, function(f){
        subset.snap.gmat.h5ad.single(
          barcodes = barcodes, h5adFile = f, colnms = genes)
      }, mc.cores = ncore)
  } else {
    lapply(sampleFiles, function(f) {
        subset.snap.gmat.h5ad.single(
          barcodes = barcodes, h5adFile = f, colnms = genes)
    })
  }
  
  left <- !sapply(matList, is.null)
  if(sum(left) < 1) {
    stop("No barcodes found in all the sample files.")
  }
  matList <- matList[left]
  mat <- do.call(rbind, matList)
  ## remote this
  ## - we have colnames and rownames, which will be kept
  ##   after rbind
  ## - there may have less barcodes
  ## rownames(mat) <- barcodes
  ## colnames(mat) <- genes
  if(length(barcodes) != nrow(mat)) {
    message(paste(nrow(mat), "over", length(barcodes), "found."))
  } else {
    message("All the barcodes are found.")
  }
  return(mat)
}

#' @export
subset.snap.gmat.h5ad.single <- function(barcodes, h5adFile, colnms = NULL) {
  if (!file.exists(h5adFile)){
    stop(h5adFile, " does not exist.")
  }
  message("working on: ", basename(h5adFile))
  mat <- h5adToMat.Sparse(filenm = h5adFile)
  index <- barcodes %in% rownames(mat)
  if(sum(index) < 1) {
    message("No barcodes found.")
    return(NULL)
  }
  r <- mat[barcodes[index], , drop = FALSE]
  if (!is.null(colnms)) {
    colnames(r) <- colnms
  }
  return(r)
}
