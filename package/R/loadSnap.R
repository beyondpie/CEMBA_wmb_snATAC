library(stringr)
## need to load cembav2env from cembav2env.R
source(file.path(here::here(), "package", "R", "cembav2env.R"))

#' @export
loadL2Snap <- function(cembav2env, class = "NonN") {
  snapL2Dir <- cembav2env$snapL2Dir
  r <- cembav2env$L1CMeta$L1Reso[
    cembav2env$L1CMeta$L1Label == class]
  snap <- readRDS(
    file.path(snapL2Dir, class,
              stringr::str_glue("{class}.k-50.km-RANN.snap.rds")))
  clusters <- read.csv(
    file = file.path(snapL2Dir, class,
                     stringr::str_glue("{class}.k-50.km-RANN.r-{r}.clusterid.txt")),
    header = FALSE, quote = "")$V1
  snap@cluster <- as.factor(clusters)
  ## cluster index from 1
  snap@metaData$cluster <- as.integer(clusters)
  snap@metaData$L2 <- paste(class, snap@metaData$cluster, sep = "_")
  snap@metaData$CellID <- paste(snap@sample, snap@barcode, sep = ".")
  return(snap)
}

#' @export
loadL3Snap <- function(cembav2env,
                       snapL2 = NULL,
                       group = "NonN_2",
                       class = "NonN") {
  snapL3Dir <- cembav2env$snapL3Dir
  if(class == "NonN"){
    l3meta <- cembav2env$NonNL3CMeta
  }
  if(class == "NonNExtra") {
    l3meta <- cembav2env$NonNL3CMetaExtra
  }
  if(class == "GLUT") {
    l3meta <- cembav2env$GLUTL3CMeta
  }
  if(class == "GABA") {
    l3meta <- cembav2env$GABAL3CMeta
  }
  r <- l3meta$L3Reso[l3meta$L2 == group]
  if (r < 0.1) {
    message("No L3 clustering for ", group)
    if (is.null(snapL2)) {
      stop("Need to provide snapL2 object for this case.")
    }
    groups <- as.vector(
      stringr::str_split(group, pattern = "-", simplify = T))
    snap <- snapL2[snapL2@metaData$L2 %in% groups, , drop = FALSE]
    clusters <- rep(1, nrow(snap))
  } else {
    snap <- readRDS(
      file.path(snapL3Dir, 
                stringr::str_glue("{group}.k-50.km-RANN.snap.rds")))
    clusters <- read.csv(
      file = file.path(snapL3Dir, 
                       stringr::str_glue("{group}.k-50.km-RANN.r-{r}.clusterid.txt")),
      header = FALSE, quote = "")$V1
  }
  snap@cluster <- as.factor(clusters)
  ## cluster index from 1
  snap@metaData$cluster <- as.integer(clusters)
  snap@metaData$L3 <- paste(group, snap@metaData$cluster, sep = "_")
  snap@metaData$CellID <- paste(snap@sample, snap@barcode, sep = ".")
  return(snap)
}

#' @export
loadL4SnapList <- function(dataDir,
                       rvec,
                       snapAll = NULL,
                       group = "NonN_2",
                       L4Col = "L4",
                       clusterCol = "L4Id") {
  snapList <- lapply(seq_along(rvec), function(i) {
    r <- rvec[i]
    g <- paste(group, i, sep = "_")
    if (r < 0.1) {
      message("No L3 clustering for ": group)
      if (is.null(snapAll)) {
        stop("Need to provide snapAll object for this case.")
      }
      groups <- as.vector(stringr::str_split(group, pattern = "-", simplify = T))
      s <- snapL2[snapAll@metaData$L2 %in% groups, , drop = FALSE]
      clusters <- rep(1, nrow(s))
    } else {
      snap <- readRDS(
        file.path(dataDir, 
                  stringr::str_glue("{g}.k-50.km-RANN.snap.rds")))
      clusters <- read.csv(
        file = file.path(dataDir, 
                         stringr::str_glue("{g}.k-50.km-RANN.r-{r}.clusterid.txt")),
        header = FALSE, quote = "")$V1
    }
    snap@cluster <- as.factor(clusters)
    ## cluster index from 1
    snap@metaData[, clusterCol] <- as.integer(clusters)
    snap@metaData[, L4Col] <- paste(g, as.integer(clusters), sep = "_")
    snap@metaData$CellID <- paste(snap@sample, snap@barcode, sep = ".")
    return(snap)
  })
  names(snapList) <- paste(group, seq_along(rvec), sep = "_")
  return(snapList)
}

#' @export
loadCEMBAv1Meta <- function(f = cembav2env$cembav1MetaFile) {
    r <- as.data.frame(data.table::fread(f, header = TRUE, sep = "\t"))
    return(r)
}

#' @export
loadL2SnapWithGmat <- function(gmatDir,
                               L2,
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
  return(snap)
}


createMergedSnap <- function(snapfullpath, ncores = 1) {
  snapPaths <- read.table(snapfullpath, header = F, quote = "")[, 1]
  snapList <- parallel::mclapply(snapPaths, function(f) {
    snm <- gsub(".snap", "", basename(f))
    message(snm, ":", f)
    snap <- SnapATAC::createSnap(file = f, sample = snm)
    return(snap)
  }, mc.cores = ncores)
  message("Merge all the snaps.")
  r <- SnapATAC::snapListRbind(snapList = snapList, checkSnap = T)
  return(r)
}
