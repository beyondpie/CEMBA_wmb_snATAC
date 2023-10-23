#' @export
getRegionPal <- function(brainRegionMeta, regionType = "MajorRegion") {
  pal <- brainRegionMeta[, paste0(regionType, "Color")]
  names(pal) <- brainRegionMeta[, regionType]
  uRegion <- unique(names(pal))
  uPal <- pal[uRegion]
  return(uPal)
}

#' @export
getRegionInfo <- function(snap,
                          brainRegionMeta,
                          groupBy = c("NonN_2"),
                          groupType = "L2") {
  smeta <- snap@metaData
  rowIndex <- smeta[, groupType] %in% groupBy
  samples <- snap@sample[rowIndex]
  barcodes <- snap@barcode[rowIndex]
  cellIds <- paste(samples, barcodes, sep = ".")
  dissects <- stringr::str_split(samples, pattern = "_",
                                 simplify = T)[,2]
  sm <- smeta[rowIndex, ]
  regionInfo <- brainRegionMeta[dissects, ]
  r <- data.frame(
    cellId = cellIds,
    dissect = dissects
  )
  r <- cbind(r, sm, regionInfo)
  rownames(r) <- r$cellId
  return(r)
}


getDissection <- function(sample) {
  d <- vapply(sample, function(s) {
    return(strsplit(s, split = "_", fixed = FALSE)[[1]][2])
  }, FUN.VALUE = "3C")
  return(d)
}

#' Get timestamp and dissect for each sample.
#' @param samples vector of characters, each element would be like
#' "CEMBA210101_1A" style.
#' @return data.frame, samples as rownames, "time" and "dissect" as
#' colnames.
#' @export
getTimeDissect <- function(samples) {
  timeDissect <- vapply(samples, function(s) {
    a <- strsplit(s, split = "_", fixed = TRUE)[[1]]
    return(c(gsub("CEMBA", "", a[1]), a[2]))
  }, FUN.VALUE = c(0, "a"))
  timeDissect <- t(timeDissect)
  if(is.null(names(samples))) {
    warning("Input has no names, so will not set rownames of output.")
  } else {
    message("Set rownames of output as names of samples.")
    rownames(timeDissect) <- names(samples)
  }
  colnames(timeDissect) <- c("time", "dissect")
  return(timeDissect)
}

#' Stat the number of both early and later replicates in each
#' dissection the meta data has.
#' @param m data.frame, meta data of time, dissect.
#' @param dm data.frame, dissection meta
#' @return numeric matrix, size [2, num_dissect]
#' @export
statReplicates <- function(m, dm) {
  udissect <- unique(m$dissect)
  nreps <- vapply(udissect, function(d) {
    times <- as.integer(with(m, time[dissect == d]))
    c(
      sum(times %in% as.character(dm[d, 1])),
      sum(times %in% as.character(dm[d, 2]))
    )
  }, FUN.VALUE = c(1, 1))
  rownames(nreps) <- c("early", "later")
  return(nreps)
}

#' @export
calJSD <- function(annot,
                   regionLevel = "MajorRegion",
                   cellLevel = "L2Annot") {
  t1 <- table(annot[, cellLevel], annot[, regionLevel])
  colSum <- colSums(t1)
  pbg <- colSum / sum(colSum)

  t2 <- t(t1)
  p <- t(scale(t2, center = FALSE, scale = colSums(t2)))

  jsd <- sapply(1:nrow(p), function(i) {
    t3 <- rbind(p[i, ], pbg)
    philentropy::JSD(x = t3, unit = "log2")
  })
  names(jsd) <- rownames(p)
  return(jsd)
}

#' @export
mouseBrainRegion <- list(
  Telencephalon = c("Isocortex", "HPF", "OLF", "AMY", "STR", "PAL"),
  Diencephalon = c("TH", "HY"),
  Midbrain = "MB",
  Hindbrain = c("Pons", "MY"),
  Cerebellum = "CB"
)

#' @export
majorRegion2Region <- data.frame(
  majorRegion = c("Isocortex", "HPF", "OLF", "AMY", "STR", "PAL",
    "TH", "HY",
    "MB",
    "Pons", "MY",
    "CB"),
  region = c(rep("Telencephalon", 6),
    rep("Diencephalon", 2),
    "Midbrain",
    rep("Hindbrain", 2),
    "Cerebellum"
  )
)
rownames(majorRegion2Region) <- majorRegion2Region$majorRegion
