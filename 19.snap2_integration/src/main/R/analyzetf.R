library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratObject)

library(stringr)
library(purrr)
# support R new feature in 4.2
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

library(R6)
library(future)
options(future.globals.maxSize = 10e9)
library(future.apply)
plan(multicore, workers = 4)

library(ggplot2)

# load my scripts
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(
  .from = "ggtheme.R",
  .directory = rdir,
  dotplotTheme)

mytheme <- dotplotTheme(legend.pos = "right")

# * class to summarize transfer label results.
TransferLabelSum <- R6::R6Class(
  classname = "TransferLabelSum",
  public = list(
    # params
    kAnchor = 5,
    feaName = "vf",
    dsNum = 50,
    allenTech = "10xv3",
    group = "neuron",
    method = "cca",
    anchor = NULL,
    tflabel = NULL,
    refLabelCols = c("subclass_id", "supertype_id", "cl"),
    tfnmprefix = function() {
      paste(self$group,
        paste0("atac.", self$allenTech),
        self$method,
        self$feaName,
        self$kAnchor,
        sep = "_")
    },
    initialize = function(kAnchor = 5,
                          feaName = "vf",
                          dsNum = "50",
                          allenTech = "10xv3",
                          group = "neuron",
                          method = "cca") {
      self$kAnchor <- kAnchor
      self$feaName <- feaName
      self$dsNum <- dsNum
      self$allenTech <- allenTech
      self$group <- group
      self$method <- method
    },
    loadAnchor = function(fromdir) {
      f <- file.path(fromdir,
        paste0(self$tfnmprefix(), ".tf.anchor.rds"))
      if(!file.exists(f)) {
        stop(f, " does not exist.")
      }
      self$anchor <- readRDS(f)
      invisible(self)
    },
    loadtfLabel = function(fromdir) {
      f <- file.path(fromdir,
        paste0(self$tfnmprefix(), ".tf.label.rds"))
      if(!file.exists(f)) {
        stop(f, " does not exist.")
      }
      self$tflabel <- readRDS(f)
      if (is.null(names(self$tflabel))) {
        message("name of tflabel is null, then add the name.")
        names(self$tflabel) <- self$refLabelCols
      }
      invisible(self)
    }
  )
)

# * functions
toNamedArray.1dtable <- function(t) {
  tmp <- as.data.frame(t, stringsAsFactors = FALSE)
  r <- tmp$Freq
  names(r) <- tmp$Var1
  return(r)
}

to3.matrix <- function(mat,
                       names = c("row", "column", "score"),
                       int2str = TRUE,
                       factor2str = TRUE) {
  r <- reshape2::melt(mat)
  colnames(r) <- names
  if (factor2str) {
    if (is.factor(r[,1])) {
      r[,1] <- as.character(r[,1])
    }
    if (is.factor(r[,2])) {
      r[,2] <- as.character(r[,2])
    }
  }
  if(int2str) {
    if(is.numeric(r[,1])) {
      r[,1] <- as.character(r[,1])
    }
    if(is.numeric(r[,2])) {
      r[,2] <- as.character(r[,2])
    }
  }
  return(r)
}

getMaxColScoreWithName <- function(mat) {
  maxScores <- apply(mat, 2, max)
  rnames <- rownames(mat)
  maxNames <- apply(mat, 2, \(x){
    rnames[which.max(x)]
  })
  names(maxScores) <- maxNames
  return(maxScores)
}

getPredictLabel <- function(tfquery) {
  labels <- tfquery$predicted.id
  scores <- tfquery$predicted.id.score
  cells <- colnames(tfquery)
  r <- data.frame(
    barcode = cells,
    predict = labels,
    score = scores
  )
  rownames(r) <- cells
  return(r)
}

getPredictScoreMat <- function(tfquery) {
  predictAssay <- tfquery@assays$prediction.score.id
  score <- predictAssay@data
  return(score)
}

#' @return array with name
getMetaCol.Seurat <- function(seu, colnm) {
  ## use repeat colnm to get array
  ## instead of data.frame with one column
  g <- seu[[colnm]][[colnm]]
  names(g) <- colnames(seu)
  return(g)
}

## TLS shor for transfer label score
getAvgTLSMat.queryGroup <- function(tfSum,
                                    refCol = "cl",
                                    queryCol = "L4",
                                    rowfn = rowSums) {
  tfmat <- getPredictScoreMat(tfSum$tflabel[[refCol]])
  nrefCol <- nrow(tfmat)
  nquery <- ncol(tfmat)
  message(nrefCol, " groups in ref col: ", refCol)
  message(nquery, " cells are predicted.")
  queryGroup <- getMetaCol.Seurat(tfSum$tflabel[[refCol]])
  ug <- unique(queryGroup)
  message("Find ", length(ug), " groups from ", queryCol)
  r <- future.apply::future_vapply(
    X = ug,
    FUN = \(g) {
      n <- sum(queryGroup == g)
      rowfn(tfmat[, queryGroup == g, drop = FALSE]) / n
    },
    FUN.VALUE = rep(0.0, nrefCol)
  )
  colnames(r) <- ug
  return(r)
}

getAvgVoteMat.queryGroup <- function(tfSum,
                                     refCol = "cl",
                                     queryCol = "L4") {
  seu <- tfSum$tflabel[[refCol]]
  tfmat <- getPredictScoreMat(seu)
  refGroup <- rownames(tfmat)
  nrefCol <- nrow(tfmat)
  nquery <- ncol(tfmat)
  message(nrefCol, " groups in ref col: ", refCol)
  message(nquery, " cells are predicted.")
  
  tfLabels <- getPredictLabel(seu)
  
  queryGroup <- getMetaCol.Seurat(seu, queryCol) |>
    x => x[tfLabels$barcode]
  
  ug <- unique(queryGroup)
  message("Find ", length(ug), " groups from ", queryCol)
  r <- future.apply::future_vapply(
    X = ug,
    FUN = \(g)  {
      r <- rep(0.0, length(refGroup))
      names(r) <- refGroup
      index <- queryGroup == g
      n <- sum(index)
      tmp <- table(tfLabels[index, "predict"]) |>
        toNamedArray.1dtable()
      r[names(tmp)] <- tmp / n
      return(r)
    },
    FUN.VALUE = rep(0.0, length(refGroup))
  )
  colnames(r) <- ug
  rownames(r) <- refGroup
  return(r)
}

#' tfmat: transfer label score: ref by query
prepareDotPlot4TransferLabel <- function(tfmat,
                                         refOrder = NULL,
                                         names = c("row", "column", "score"),
                                         ignoreEmptyRef = TRUE) {
  maxScore <- getMaxColScoreWithName(mat = tfmat)
  query2ref <- data.frame(
    query = colnames(tfmat),
    ref = names(maxScore),
    row.names = colnames(tfmat)
  )
  if (ignoreEmptyRef) {
    message("remove refs not having query mapped to.")
    tfmat <- tfmat[rownames(tfmat) %in% query2ref$ref, ]
  }
  if (is.null(refOrder)) {
    message("refOrder is null, will use default numeric order for it.")
    refOrder <- rownames(tfmat)
    refOrder <- refOrder[order(as.integer(refOrder))]
  } else {
    refOrder <- refOrder[refOrder %in% rownames(tfmat)]
  }
  queryOrder <- query2ref$query[
    order(factor(query2ref$ref, levels = refOrder))]

  meltmat <- to3.matrix(tfmat, names)
  meltmat[,1] <- factor(meltmat[,1], levels = refOrder)
  meltmat[,2] <- factor(meltmat[,2], levels = queryOrder)
  # reduce size of meltmat
  meltmat <- meltmat[meltmat[,3] > 0, ]
  return(meltmat)
}


plotf.hardcode <- function(fromdir,
                           outfigdir,
                           kAnchor = 5,
                           feaName = "vf",
                           dsNum = 50,
                           allenTech = "10xv3",
                           group = "neuron",
                           method = "rpca",
                           refCols = c("subclass_id", "supertype_id", "cl"),
                           queryCol = "L4",
                           ignoreEmptyRef = TRUE) {
  message("load tfSum.")
  tfSum <- TransferLabelSum$new(kAnchor = kAnchor,
    feaName = feaName, dsNum = dsNum,
    allenTech = allenTech,
    group = group,
    method = method)

  ## tfSum$loadAnchor(fromdir = fromdir)
  message("load tf label result.")
  tfSum$loadtfLabel(fromdir = fromdir)

  ## tfmat <- getAvgTLSMat.queryGroup(tfSum,
  ##   refCol = refCol,
  ##   queryCol = queryCol)

  ## tfbymax <- getMaxColScoreWithName(tfmat)
  for (refCol in refCols) {
    message("get tfVoteMat for refCol: ", refCol)
    tfVoteMat <- getAvgVoteMat.queryGroup(tfSum,
      refCol = refCol,
      queryCol = queryCol)

    message("prepare dotplot mat.")
    mat.dotplot <- prepareDotPlot4TransferLabel(
      tfmat = tfVoteMat, ignoreEmptyRef = ignoreEmptyRef,
      names = c("Allen", "ATAC", "score"))

    nref <- length(unique(mat.dotplot[, 1]))
    ntref <- nrow(tfVoteMat)
    nquery <- ncol(tfVoteMat)

    lowSimScore <- quantile(mat.dotplot$score, 0.1)
    highSimScore <- max(mat.dotplot$score)

    p.tfvote <- ggplot(data = mat.dotplot, aes(x = ATAC, y = Allen)) +
      geom_point(aes(color = score, size = score)) +
      mytheme +
      theme(
        panel.background = element_rect(fill = "#132B43"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
      scale_color_gradient(low = "#132B43", high = "#56B1F7",
        limits = c(lowSimScore, highSimScore), na.value = "#132B43") +
      scale_size(range = c(0, 1)) +
      xlab(str_glue("{nquery} {queryCol} of ATAC")) +
      ylab(str_glue("{nref}/{ntref} {refCol} of Allen")) +
      ggtitle(tfSum$tfnmprefix())

    message("save figure.")
    ggsave(plot = p.tfvote,
      filename = file.path(outfigdir,
        paste(tfSum$tfnmprefix(), refCol, "pdf", sep = ".")),
      width = 13, height = 10)
  }
  message("Done.")
}

# * main
args <- commandArgs(trailingOnly = TRUE)
kAnchor <- as.integer(args[1])
feaName <- args[2]
allenTech <- args[3]
method <- args[4]
group <- args[5]

fromdir <- file.path(here::here(), "19.snap2_integration",
  "out/tf_tscc", group)
outfigdir <- file.path(here::here(),
  "19.snap2_integration", "out", "figure_tf")
if (!dir.exists(outfigdir)) {
  dir.create(outfigdir, recursive = TRUE)
}

# test
## plotf.hardcode(fromdir = fromdir,
##   kAnchor = 5,
##   feaName = "vf",
##   dsNum = 50,
##   allenTech = "10xv3",
##   group = "neuron",
##   method = "rpca",
##   queryCol = "L4",
##   ignoreEmptyRef = TRUE)
plotf.hardcode(fromdir = fromdir,
  outfigdir = outfigdir,
  kAnchor = kAnchor,
  feaName = feaName,
  dsNum = 50,
  allenTech = allenTech,
  group = group,
  method = method,
  queryCol = "L4",
  ignoreEmptyRef = TRUE)
