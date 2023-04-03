#' @export
colFn <- function(col.fn = min) {
  a <- function(mat) {
    apply(mat, 2, col.fn)
  }
  return(a)
}

#' @export
rowFn <- function(row.fn = min) {
  a <- function(mat) {
    apply(mat, 1, row.fn)
  }
  return(a)
}

#' Prepare path for a new output file.
#' @param outf characters, file name with desired path
#' @param remove bool, if remove old file, default TRUE
#' @return None
#' @export
cleanOutfile <- function(outf, remove = TRUE) {
  outdir <- dirname(outf)
  if (!dir.exists(outdir)) {
    message(paste(outdir, "does not exist and will be created."))
    dir.create(outdir, recursive = T)
  }
  if (file.exists(outf) & remove) {
    message(paste(outf, "exists and will be removed."))
    file.remove(outf)
  }
}

#' If file does not exists, then prepare it.
#' @param f characters, file name
#' @return None
#' @export
prepareOutfile <- function(f, remove = TRUE) {
  if (!file.exists(f)) {
    outdir <- dirname(outf)
    if (!dir.exists(outdir)) {
      message(paste(outdir, "does not exist and will be created."))
      dir.create(outdir, recursive = T)
    }
  }
  message(f, " is already existed.")
  if (remove) {
    messge("Remove ", f)
    file.remove(f)
  }
}

#' If dir does not exists, then create it.
#' @param d characters, dir name
#' @return None
#' @export
prepareOutdir <- function(d) {
  if (!dir.exists(d)) {
    message(paste(d, "does not exist and will be created."))
    dir.create(d, recursive = T)
  }
}

#' Check if element in argument has no defination
#' argument generated from package optparse
#' But this function can be used for a list
#' @param args list-like object with attribute names
#' @return None, if no defined attribute in args, stop will be called
#' @export
checkArgsExistOrStop <- function(args) {
  invisible(lapply(names(args), function(v) {
    if (is.null(args[[v]])) {
      stop("Args have no attribute ", v)
    }
  }))
}

#' Check if file/dir exist, if not, then stop
#' @param f characters, a file or director
#' @return None, will send stop signal if
#' f doesn't exist.
#' @export
checkFileExistOrStop <- function(f) {
  if ((!file.exists(f)) & (!dir.exists(f))) {
    stop(f, " does not exist.")
  }
}

#' @export
split1 <- function(str, split = "-") {
  return(unlist(strsplit(str, split = split, fixed = TRUE)))
}

#' @export
split2 <- function(str, s1 = "-", s2 = "_", pos = NULL){
  r1 <- split1(str, split = s1)
  r2 <- lapply(r1, split1, split = s2)
  names(r2) <- paste("s", seq_len(length(r2)), sep = "")
  if(is.null(pos)) {
    return(r2)
  }
  return(vapply(r2, function(i){i[pos]}, FUN.VALUE = "c1"))
}

#' @export
getcov.array <- function(x) {
  cv <- sd(x) / abs(mean(x))
  return(cv)
}


#' Get coefficicent of variation for each features.
#' mat: cell by features
#' @export
getcov.mat <- function(mat) {
  cov <- apply(mat, 2, getcov.array)
  return(cov)
}

#' @export
setNamedVector <- function(df, kcol, vcol) {
  tmp <- unique(df[, c(kcol, vcol)])
  v <- tmp[, vcol]
  names(v) <- tmp[, kcol]
  return(v)
}

#' @export
sparse2triples <- function(m) {
 SM = summary(m)
 D1 = m@Dimnames[[1]][SM[,1]]
 D2 = m@Dimnames[[2]][SM[,2]]
 data.frame(peak=D2, cell=D1, x=m@x)
}

#' @export
shuf.mat.col <- function(mat, seed = 2023) {
  set.seed(seed)
  r <- mat[ , sample(ncol(mat))]
  colnames(r) <- colnames(mat)
  rownames(r) <- rownames(mat)
  return(r)
}

#' @export
shuf.mat.row <- function(mat, seed = 2023) {
  set.seed(seed)
  r <- mat[sample(nrow(mat)), ]
  colnames(r) <- colnames(mat)
  rownames(r) <- rownames(mat)
  return(r)
}

#' @export
aggMatByCol <- function(mat, group.col, fn = rowSums) {
  ug <- unique(group.col)
  r <- vapply(ug, function(g) {
    t <- fn(mat[, group.col %in% g])
  }, FUN.VALUE = rep(0.0, nrow(mat)))
  colnames(r) <- ug
  return(r)
}


#' @export
scaleByRow <- function(mat,
                       is.center = TRUE,
                       is.scale = TRUE,
                       rmZeroSdRows = TRUE,
                       is.capped = TRUE,
                       up.hard = NULL,
                       down.hard = NULL,
                       up.quantile = 0.999,
                       down.quantile = 0.001
                       ) {
  sdRows <- apply(mat, 1, sd)
  zeroRows <- (sdRows ==0)
  if(any(zeroRows)) {
    message(sum(zeroRows), " has zero standard deviation.")
    if (rmZeroSdRows) {
      message("Will remove them.")
      mat <- mat[!zeroRows, ]
    } else {
      message("Data may have nan or NA after scaling.")
    }
  }
  r <- t(scale(t(mat), center = is.center, scale = is.scale))
  rownames(r) <- rownames(mat)
  colnames(r) <- colnames(mat)
  if (is.capped) {
    message("Will cap the data")
    up.thres <- if(!is.null(up.hard)) {
      up.hard
    } else {
      quantile(r, up.quantile)
    }
    down.thres <- if(!is.null(down.hard)) {
      down.hard
    } else {
      quantile(r, down.quantile)
    }
    message(
      paste("Restrict data between", down.thres, "and", up.thres, sep = " "))
    r[r >= up.thres] <- up.thres
    r[r <= down.thres] <- down.thres
  }
  return(r)
}


#' @export
scaleUseRankByRow <- function(mat, tiesMethod = "average") {
  rank_mat <- apply(mat, 1, rank, ties.method = tiesMethod)
  r <- scaleByRow(mat = t(rank_mat), is.center = TRUE,
    is.scale = TRUE, rmZeroSdRows = TRUE,
    is.capped = FALSE)
  return(r)
}

#' @param mat1, row scaled mat, same shape as mat2
#' @param mat2, row scaled mat, same shape as mat1
fastPearCorOnRow.scaled <- function(mat1, mat2) {
  r <- rowSums(mat1 * mat2) / ncol(mat1)
  # ideally we do not need this
  # but our mat may be use hard threshold to
  # limit the too high / too low values
  # So add this just in case
  r[r >=1 ] <- 1
  r[r <= (-1)] <- (-1)
  return(r)
}

#' @export
#' @examples
#'
#' n <- 20000
#' row.index <- sample(min(nrow(s.rna.mat),nrow(s.atac.mat)), n)
#' ## pearson 
#' r1 <- fastCorrOnRow(mat1 = s.rna.mat[row.index, ], mat2 = s.atac.mat[row.index, ])
#' r2 <- vapply(row.index, function(i) {
#'   cor(x = s.rna.mat[i,], y = s.atac.mat[i,], method = "pearson")
#' }, FUN.VALUE = 0.1)
#' plot(x = r1, y = r2, xlab = "factCor", ylab = "cor")
#' abline(a=0, b=1)
#' qqplot(x = r1, y = r2)
#' abline(a=0, b=1)
#' ## spearman
#' r1 <- fastCorrOnRow(
#'   mat1 = normalizeUseRankByRow(s.rna.mat[row.index, ]),
#'   mat2 = normalizeUseRankByRow(s.atac.mat[row.index, ]))
#' r2 <- vapply(row.index, function(i) {
#'   cor(x = s.rna.mat[i,], y = s.atac.mat[i,], method = "spearman")
#' }, FUN.VALUE = 0.1)
#' plot(x = r1, y = r2, xlab = "factCor", ylab = "cor")
#' abline(a=0, b=1)
#' qqplot(x = r1, y = r2)
#' abline(a=0, b=1)
#' 
#' @param mat1 same shape as mat2
#' @param mat2 same shape as mat1
fastSpearCorOnRow.distinctRank <- function(mat1, mat2, tiesMethod = "average") {
  rankMat1 <- t(apply(mat1, 1, rank, ties.method = tiesMethod))
  rankMat2 <- t(apply(mat2, 1, rank, ties.method = tiesMethod))
  n <- ncol(mat1)
  diffRank <- 1 -  6 * rowSums((rankMat1-rankMat2)^2) / (n * (n^2-1))
  return(diffRank)
}

#' @param p1byc scaled mat with row and col names
#' ncol(p1byc) = ncol(p2byc)
#' @param p2byc scaled mat with row and col names
#' ncol(p1byc) = ncol(p2byc)
#' @export
fastCorOnRow.chunk <- function(p1byc, p2byc, p12pairs,
                                chunkSize = 10e+05,
                                corMethod = c("spearman", "pearson")) {
  method <- match.arg(corMethod)
  # * filter p12pairs based on rownames of p1byc, p2byc
  upairs <- unique(p12pairs)
  message("Unique pairs: ", nrow(upairs))
  p1s <- unique(upairs[,1])
  p2s <- unique(upairs[,2])
  p1s.used <- p1s[p1s %in% rownames(p1byc)]
  p2s.used <- p2s[p2s %in% rownames(p2byc)]
  message(paste(length(p1s.used), "/", length(p1s),"ids are in p1byc."))
  message(paste(length(p2s.used), "/", length(p2s),"ids are in p2byc."))
  
  uupairs <- upairs[(upairs[,1] %in% p1s.used) & (upairs[,2] %in% p2s.used), ]
  message("After filtering names: ", nrow(uupairs), " pairs left.")

  # * use common columns in p1byc, p2byc
  com.cols <- intersect(colnames(p1byc), colnames(p2byc))
  message(length(com.cols), " common columns in p1byc and p2byc.")
  p1byc <- p1byc[, com.cols]
  p2byc <- p2byc[, com.cols]
  
  index.all <- seq_len(nrow(uupairs))
  chunks <- split(index.all, ceiling(seq_along(index.all) / chunkSize))
  if (method == "pearson") {
    message("Remind: p1byc and p2byc should be scaled for pearson.")
  }
  if (method == "spearman") {
    message("Remind: p1byc and p2byc should be scaled by rank for spearman")
    ## p1byc <- scaleUseRankByRow(p1byc)
    ## p2byc <- scaleUseRankByRow(p2byc)
  }
  rs <- lapply(chunks, function(index.chunk) {
    p12 <- uupairs[index.chunk, ]
    mat1 <- p1byc[p12[,1], ]
    mat2 <- p2byc[p12[,2], ]
    fastPearCorOnRow.scaled(mat1 = mat1, mat2 = mat2)
  })
  message("Merge results from all the chunks.")
  r <- unlist(rs)
  names(r) <- paste(uupairs[,1], uupairs[,2], sep = "@")
  return(r)
}

#' @export
df2tri <- function(df, byrow = TRUE) {
  if(!byrow) {
    df <- t(df)
  }
  triList <- lapply(rownames(df), function(i) {
    data.frame(row = rep(i, ncol(df)),
      col = colnames(df), val = unlist(df[i, ]))
  })
  do.call(rbind, triList)
}
