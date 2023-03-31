# library(S4Vectors)
#' @export
initAnnot <- function(cemballMeta) {
  cemballMeta$L1Annot <- cemballMeta$L1
  cemballMeta$L2Annot <- cemballMeta$L2
  cemballMeta$L3Annot <- cemballMeta$L3
  return(cemballMeta)
}

#' @param allmeta meta from snap@metaData.
#' @param l1 characters, like "NonN", "GLUT" or "GABA"
#' @param l3meta is from cemba2env GLUT or NonN L3CMeta
#' @export
#' @importFrom stringr str_glue
#' @importFrom S4Vectors first second
addL3Annot <- function(allmeta,
                       l1,
                       l3meta) {
  result <- allmeta
  L2 <- l3meta$L2
  L2Annot <- l3meta$L2Label
  names(L2Annot) <- L2
  L3n <- l3meta$L3nCluster
  names(L3n) <- L2
  L3Annot <- l3meta$L3Label
  resetL1 <- l3meta$resetL1
  resetL2 <- l3meta$resetL2
  result$L1Annot[result$L1 %in% l1] <- l1
  # * set l2annot
  for (l2elem in L2) {
    result$L2Annot[result$L2 %in% l2elem] <-
      L2Annot[l2elem]
  }
  # * set l3annot
  notInL3Annot <- !(L2 %in% names(L3Annot))
  if (sum(notInL3Annot) > 0) {
    message(stringr::str_glue(
      "{sum(notInL3Annot)} L2 does not have L3 annot."))
    message("They are: ", paste(L2[notInL3Annot], collapse = ","))
  } else {
    message("All the L2s have L3 annot record.")
  }
  for (l2elem in names(L3Annot)) {
    message(stringr::str_glue("set L3 annot for L2: {l2elem}."))
    l3annot <- L3Annot[[l2elem]]
    if (any(duplicated(S4Vectors::first(l3annot)))) {
      stop("No unique element in L3 annot.")
    }
    nl3a <- length(l3annot)
    nl3 <- L3n[l2elem]
    if (nl3a <= nl3) {
      message(stringr::str_glue("{nl3a} of {nl3} have L3 annot."))
    } else {
      stop(stringr::str_glue("{nl3a} L3 annot > {nl3} clusters we have."))
    }
    for (i in seq_along(l3annot)) {
      i1 <- paste(l2elem, S4Vectors::first(l3annot[i]), sep = "_")
      i2 <- S4Vectors::second(l3annot[i])
      message(stringr::str_glue("Annot {i1} as {i2}."))
      result$L3Annot[result$L3 %in% i1] <- i2
    }
  }
  # ** reset l2annot
  if (!is.null(resetL2)) {
    message("resetL2 is not null. update L2Annot.")
    for (i in seq_along(resetL2)) {
      i1 <- S4Vectors::first(resetL2[i])
      i2 <- S4Vectors::second(resetL2[i])
      message(stringr::str_glue("Update {i1}'s class as {i2}."))
      result$L2Annot[grepl(i1, result$L3Annot)] <- i2
    }
  }
  # ** reset l1annot
  if (!is.null(resetL1)) {
    message("resetL1 is not null. update L1Annot.")
    for (i in seq_along(resetL1)) {
      i1 <- S4Vectors::first(resetL1[i])
      i2 <- S4Vectors::second(resetL1[i])
      message(stringr::str_glue("Update {i1}'s class as {i2}."))
      result$L1Annot[grepl(i1, result$L2Annot)] <- i2
    }
  }
  return(result)
}

# due to time, we assume all the subclusters have annotations.
#' @export
updateAnnotWithL3Extra <- function(cemballMeta,
                                   l1 = "NonN",
                                   l3meta) {
  L2 <- l3meta$L2
  L3n <- l3meta$L3nCluster
  names(L3n) <- L2
  L3Label <- l3meta$L3Label
  L4Label <- l3meta$L4Label
  for (i in seq_along(L2)) {
    l2elem <- L2[i]
    l2ids <- gsub(paste0(l1, "_"), "",
      as.vector(stringr::str_split(l2elem, pattern = "-",
        simplify = TRUE)))
    l3nm <- paste(l1, paste(l2ids, collapse = "."), sep = "_")
    l3Label <- L3Label[[l3nm]]
    for (j in seq_along(l3Label)) {
      i1 <- S4Vectors::first(l3Label[j])
      i2 <- S4Vectors::second(l3Label[j])
      l3extra <- paste(l2elem, i1, sep = "_")
      cemballMeta$L2Annot[
        cemballMeta$L3Extra %in% l3extra] <- i2
      l4nm <- paste(l3nm, i1, sep = "_")
      if (l4nm %in% names(L4Label)) {
        l4Label <- L4Label[[l4nm]]
        for (t in seq_along(l4Label)) {
          i1 <- S4Vectors::first(l4Label[t])
          i2 <- S4Vectors::second(l4Label[t])
          cemballMeta$L3Annot[
            cemballMeta$L4 %in% stringr::str_glue("{l3extra}_{i1}")] <- i2
        } # end of t
      } # endf of if l4nm
    } # end of j
  } # end of i
  return(cemballMeta)
}

#' @export
setClassForSubclass <- function(meta = NULL,
                                cembav2env = NULL,
                                classCol = "upL1Annot",
                                subclassCol = "upL2Annot") {
  if (is.null(meta)) {
    if (is.null(cembav2env)) {
      source(file.path(here::here(), "package", "R", "cembav2env.R"))
    }
    meta <- readRDS(cembav2env$allmetaFile)
  }
  class2subclass <- unique(meta[, c(classCol, subclassCol)])
  colnames(class2subclass) <- c("class", "subclass")
  class2subclass$class[
    class2subclass$subclass %in% "EPEN"] <- "NN"
  class2subclass$class[
    class2subclass$subclass %in% "Hmx2DMHGA"] <- "Gaba"
  class2subclass$class[
    class2subclass$subclass %in% "SVRDOPA"] <- "Dopa"
  class2subclass <- unique(class2subclass)
  return(class2subclass)
}

#' @export
setUClassForSubclass <- function(class) {
  uclass <- class
  uclass[grep("NN|nn", uclass)] <- "NN"
  uclass[grep("Gaba;Glut", uclass)] <- "Glut"
  uclass[grep("Chol", uclass)] <- "Chol"
  uclass[grep("DOPA|Dopa", uclass)] <- "Dopa"
  uclass[grep("Gly", uclass)] <- "Gly"
  uclass[grep("Nora", uclass)] <- "Nora"
  uclass[grep("Glut;Gaba", uclass)] <- "Gaba"
  return(uclass)
}

#' Order subclasses based on class information.
#' NN, Glut, Gaba mainly
#' @export
order.subclass.1 <- function(subclass, class) {
  class <- toupper(class)
  subclass.nn <- subclass[grepl("NN", class)]
  subclass.glut <- subclass[grepl("GLUT", class)]
  ## So Gaba;Glut or Glut;Gaba will be labeled as GLUT.
  subclass.gaba <- setdiff(
    x = subclass, y = c(subclass.nn, subclass.glut))
  r <- data.frame(
    subclass = c(subclass.nn, subclass.glut, subclass.gaba),
    mainclass = c(
      rep("NN", length(subclass.nn)),
      rep("GLUT", length(subclass.glut)),
      rep("GABA", length(subclass.gaba))
  ))
  return(r)
}


#' @export
getSubclassMeta.ordered <- function(file,
                                    col.subclass = "subclass.hc",
                                    orderby = "index.hc.byclass",
                                    col.class = "uclass",
                                    color.col.class = "uclassColor",
                                    color.col.subclass = "subclasscolor",
                                    which.ha = "row",
                                    legend.class = TRUE,
                                    legend.subclass = FALSE) {
  subclassMeta <- data.table::fread(
    file = file, header = TRUE, sep = ",", data.table = FALSE)
  subclass.nm <- subclassMeta[, col.subclass]
  rownames(subclassMeta) <- subclass.nm
  r <- subclassMeta[subclassMeta[, orderby], ]
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Class = r[, col.class],
    Subclass = r[, col.subclass],
    which = which.ha,
    col = list(
      Class = setNames(unique(r[, color.col.class]),
        unique(r[, col.class])),
      Subclass = setNames(unique(r[, color.col.subclass]),
        unique(r[, col.subclass]))
    ),
    show_legend = c(legend.class, legend.subclass)
  )
  return(list(r = r, ha = ha))
}


#' @export
concatAllenAnnot <- function(allenAnnot,
                             space2 = "_",
                             slash2 = "-") {
  r <- gsub("/", slash2, allenAnnot)
  r <- gsub(" ", space2, r)
  return(r)
}

#' Deprecated
#' @export
mapAllenSubclassTomCGSubclass <- function(allenClass) {
  r <- gsub(pattern = " ", replacement = "_", allenClass)
  r1 <- gsub("/", "_", r)
  r2 <- gsub("_$", "", r1)
  return(r2)
}

#' @export
mapAllenSubclassTomCGSubclass2 <- function(allenClass) {
  r <- gsub(pattern = " ", replacement = "_", allenClass)
  r1 <- gsub("/", "-", r)
  r2 <- gsub("_$", "", r1)
  return(r2)
}
