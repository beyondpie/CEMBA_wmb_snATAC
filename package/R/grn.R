#' @export
loadRawGRNFile <- function(f, indir = NULL) {
  message("Loading ", f)
  if(!is.null(indir)) {
    f <- file.path(indir, f)
  }
  d <- data.table::fread(file = f, header = TRUE, sep = ",",
    data.table = FALSE)
  r <- d[ , c(2,3,4,7)]
  colnames(r) <- c("source", "target", "coef", "nlogp")
  return(r)
}

#' @export
getSubclassFromGRNFiles <- function(GRNFiles) {
  subclasses <- vapply(GRNFiles, function(i) {
    unlist(strsplit(i, split = ".", fixed = TRUE))[2]
  }, FUN.VALUE = "ZI_Pax6_Gaba")
  return(subclasses)
}

#' @export
filterGRN <- function(rawGRNs, lnlogp = 2, topk = 5000) {
  GRNs <- lapply(rawGRNs, function(i) {
    absCut <- sort(abs(i$coef), decreasing = TRUE)[
      min(nrow(i), topk)]
    r <- subset(i, (nlogp >= lnlogp) & (abs(coef) >= absCut) )
    return(r)
  })
}

#' @export
getAttrFromGRNs <- function(GRNs, attri = "source") {
  unique(unlist(lapply(GRNs, function(i) {i[ , attri]})))
}

#' @export
getGRNmats <- function(GRNs, TFs, Genes) {
  GRNmats <- lapply(GRNs, function(grn) {
    i <- match(grn$source, TFs)
    j <- match(grn$target, Genes)
    x <- grn$coef
    sparseMatrix(i = i, j = j, x = x,
      dims = c(length(TFs), length(Genes)),
      dimnames = list(TFs, Genes), index1 = TRUE)
  })  
}

#' Deprecated
#' Use parallel edges to represent neg loop
#' @export
toiGraph1 <- function(TF2) {
  edgeList <- lapply(seq_len(nrow(TF2)), function(i) {
    row <- TF2[i, ]
    edges <- if(row[3] < 0) {
      rep(c(row[1], row[2]), 2)
    } else {
      c(row[1], row[2])
    }
  })
  edges <- unlist(edgeList)
  make_graph(edges = edges, directed = TRUE)
}

#' Use color to label inhibit or activate network.
#' @export
toiGraph2 <- function(TF2) {
  tf2 <- TF2[, c("source", "target")]
  edges <- as.vector(t(tf2))
  g <- igraph::make_graph(edges = edges)
  igraph::E(g)$fn <- rep("activate", nrow(tf2))
  igraph::E(g)$fn[TF2$coef < 0] <- "block"
  igraph::E(g)$color <- ifelse(
    igraph::E(g)$fn == "activate", 1, 2)
  return(g)
}

#' Use colors to label none, inhibit, activate cons
#' @export
toiGraph3 <- function(TF2, none = 3, inh = 2, act = 1) {
  tf2 <- TF2
  rownames(tf2) <- with(tf2, paste(source, target, sep = ":"))
  tf2$color <- rep(act, nrow(tf2))
  tf2$color[tf2$coef < 0] <- inh
  allNodes <- unique(c(tf2$source, tf2$target))
  completeTF2 <- expand.grid(allNodes, allNodes,
    stringsAsFactors = FALSE)
  colnames(completeTF2) <- c("source", "target")
  rownames(completeTF2) <- with(completeTF2,
    paste(source, target, sep = ":"))
  completeTF2$color <- rep(none, nrow(completeTF2))
  rows.in.tf2 <- rownames(completeTF2)[
    rownames(completeTF2) %in% rownames(tf2)]
  completeTF2[rows.in.tf2, "color"] <- tf2[rows.in.tf2, "color"]
  edges <- as.vector(t(completeTF2[, 1:2]))
  g <- igraph::make_graph(edges = edges)
  igraph::E(g)$color <- completeTF2$color
  return(g)
}

#' @export
toTF2 <- function(GRN, TFs){
  r <- subset(GRN, target %in% TFs)
  return(r)
}

#' @export
countMotifs <- function(netMotifList, g) {
  vapply(netMotifList, function(m) {
    igraph::count_subgraph_isomorphisms(
      pattern = m,
      target = g,
      method = "vf2"
    )
  }, FUN.VALUE = 1)
}

#' @export
getMotifs <- function(netMotifList, g) {
  motifs <- lapply(seq_along(netMotifList), function(i) {
    m <- netMotifList[[i]]
    edgeList <- igraph::subgraph_isomorphisms(
      pattern = m,
      target = g,
      method = "vf2"
    )
    r <- if(length(edgeList) < 1){
      message("No motif found in ", names(netMotifList)[i])
      NULL
    } else {
      nmList <- lapply(edgeList, function(e) {
        names(e)
      })
      do.call(rbind, nmList)
    }
    return(r)
  })
  names(motifs) <- names(netMotifList)
  return(motifs)
}

#' Depend on toiGraph2 or toiGraph3
#' @export
sparseMatrix2igraph <- function(s, graphfn) {
  ## summary(s)
  ## or idx <- which(A != 0, arr.ind=TRUE)
  ss <- as(s, "TsparseMatrix")
  source <- ss@i + 1
  target <- ss@j + 1
  coef <- ss@x
  source_nm <- rownames(s)[source]
  target_nm <- colnames(s)[target]
  TF2 <- data.frame(source = source_nm, target = target_nm, coef = coef)
  g <- graphfn(TF2)
}

#' Depend on sparsematirx2igraph and toiGraph2
#' @export
fitPowerLaw <- function(aGRN, graphfn = toiGraph2) { 
  origG <- sparseMatrix2igraph(s = aGRN, graphfn = graphfn)
  d <- igraph::degree(origG, mode = "all")
  ## dd <- degree_distribution(graph = origG, cumulative = FALSE, mode = "all")
  fit <- igraph::fit_power_law(d)
  message(paste(
    round(c(fit$alpha, fit$KS.p, fit$xmin),
      digits = 2), collapse = ";"))
  if (fit$KS.p <= 0.05) {
    message("Power fitting failed.")
  }
  return(fit)
}

