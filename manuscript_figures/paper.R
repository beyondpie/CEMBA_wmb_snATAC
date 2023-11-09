library(patchwork)
library(tidyverse)
library(ggvenn)
library(data.table)
library(ggforce)
library(fastmatrix)
library(ComplexHeatmap)
library(Matrix)
library(future)
options(future.globals.maxSize = 5e9)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

projdir <- here::here()
packdir <- file.path(projdir, "package/R")

import::from(.from = "cembav2env.R", .directory = packdir,
  cembav2env, Sa2PeakCalling, Sa2Integration)

import::from(.from = "utils.R", .directory = packdir,
  setNamedVector, split1, split2,
  getcov.array, getcov.mat, colFn,
  scaleByRow, df2tri)

import::from(.from = "peak.R", .directory = packdir,
  loadParsePmat, processCPM.log2.scale.cap,
  orderModules.cpm.statPeak, orderModules.cpm.statPeak2,
  orderStatPeakByModuleRank, getPeakModuleAvgScore.NMF,
  getTopRankClusterForPeakModule, downsample.index.1,
  downsampleStatPeak.1, rankModules, loadStatPeak.NMF)

import::from(.from = "hc.R", .directory = packdir,
  plot_hc, gethc.classlevel, getNMFHeatMap)

import::from(.from = "integration.R", .directory = packdir,
  getIntOvlpMat)

import::from(.from = "plot.R", .directory = packdir,
  gen_circle, drawVenn3, prepareVennData.d3)

import::from(.from = "annot.R", .directory = packdir,
  getSubclassMeta.ordered)

import::from(.from = "colors.R", .directory = packdir,
  ArchRPalettes, largeRegionColors)

import::from(
  .from = "cicero.R", .directory = packdir, loadCor)

import::from(.from = "grn.R", .directory = packdir,
  toiGraph3, toiGraph2, sparseMatrix2igraph, fitPowerLaw,
  getAttrFromGRNs, toTF2, countMotifs, getMotifs)


# * functions
# ** NMF
getGlobalModule <- function(signal.scbymod,
                            quantileHighSignal =  0.95,
                            ratioOfSubclass = 0.9) {
  nsc <- nrow(signal.scbymod) * ratioOfSubclass
  highSignal <- quantile(signal.scbymod, quantileHighSignal)
  index <- colSums(signal.scbymod >= highSignal) >= nsc
  if (sum(index) < 1) {
    message("No global modules found.")
    return(NA)
  }
  return(colnames(signal.scbymod)[index])
}

setHeatmap <- function(mat,
                       low.quantile = 0.01,
                       show.legend = TRUE,
                       legend.direction = "horizontal") {
  highval <- quantile(mat, 1 - low.quantile)
  lowval <- quantile(mat, low.quantile)
  colfn <- circlize::colorRamp2(seq(lowval, highval, length = 60),
    viridis::viridis(60))
  ComplexHeatmap::Heatmap(
    matrix = mat,
    col = colfn,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    top_annotation = NULL,
    left_annotation = NULL,
    show_heatmap_legend = show.legend,
    heatmap_legend_param = list(
      title = latex2exp::TeX(r"($\log$(CPM+1))"),
      at = c(lowval, highval),
      labels = c(round(lowval, 1), round(highval, 1)),
      direction = legend.direction)
  )
}

combineHeatmap <- function(logCPM, mCG, low.quantile = 0.01) {
  hmATAC <- setHeatmap(mat = logCPM)
  hmCG <- setHeatmap(mat = mCG)
  ComplexHeatmap::draw(
    hmATAC + hmCG,
    heatmap_legend_side = "bottom",
    merge_legend = TRUE)
}

#' @param modules ordered
getMotifMat <- function(modules, motifs, motifList,
                        lowNegLogp = 10, lowAs=0,
                        globalModules = NULL,
                        highEnrichNegLogp = 50,
                        lowHighEnrichedModules = 5,
                        moduleRepeats = NULL) {
  motifMat <- matrix(data = 0, nrow = length(modules),
    ncol = length(motifs))
  rownames(motifMat) <-  modules
  colnames(motifMat) <- motifs
  for (m in modules) {
    t1 <- motifList[[m]]
    nlogp <- 0 - t1$logp
    nlogp[nlogp < lowNegLogp ] <- lowAs
    motifMat[m, t1$motif] <- nlogp
  }
  
  motif.col.order <- if(is.null(globalModules)) {
    order(motifMat[1,], decreasing = TRUE)
  } else {
    order(
      colSums(motifMat[globalModules, , drop = FALSE]),
      decreasing = TRUE, na.last = TRUE)
  }
  
  motifMat <- motifMat[, motif.col.order]
  # filter columns
  message(paste("Before filtering low motifs: ", ncol(motifMat), "motifs."))
  a <- (motifMat >= highEnrichNegLogp)
  motifMat <- motifMat[, colSums(a) >= lowHighEnrichedModules]
  message(paste("After filtering low motifs: ", ncol(motifMat), "motifs."))
  
  # expand motifmat
  r <- if (!is.null(moduleRepeats)){
    message("Expand module based on the repeats.")
    expandModules <- unlist(lapply(rownames(motifMat), function(i) {
      rep(i, moduleRepeats[i])
    }))
    motifMat[expandModules, ]
  } else {
    motifMat
  }
  return(r)  
}

getMotifModuleMat <- function(motifMat, k = 54, capped = 100) {
  d <- factoextra::get_dist(t(motifMat), method = "spearman")
  hc <- hclust(d, method = "ward.D2")
  g <- cutree(hc, k = k)
  ngroup <- length(unique(g))
  enrich <- vapply(seq_len(ngroup), function(i) {
    tmp <- motifMat[, names(g)[g == i], drop=FALSE]
    if(!is.null(capped)) {
      tmp[tmp >= capped] <- capped
    }
    r <- rowMeans(motifMat[, names(g)[g == i], drop=FALSE])
    return(r)
  }, FUN.VALUE = rep(1, nrow(motifMat)))
  colnames(enrich) <- paste0("g", seq_len(ngroup))
  return(list(enrich = enrich, group = g))
}

getOrderOfMotifMat <- function(motifMat, k = 54, topk = 1,
                               capped = 100, enrichBound = 90,
                               topGlobal = 5) {
  r <- getMotifModuleMat(motifMat, k = k, capped = capped)
  enrich <- r$enrich
  g <- r$group
  topRank.module <- getTopRankClusterForPeakModule(
    cluster2moduleScore = enrich,
    modules = colnames(enrich),
    topk = topk)
  colMins <- colFn(col.fn = min)
  modules.order <- rankModules(
    topRank.module,
    modules = colnames(topRank.module),
    avg.fn = colMins
  )
  ## nHighModule <- colSums(enrich >=enrichBound)
  ## globalModules <- names(nHighModule)[order(nHighModule, decreasing = TRUE)]
  ## globalModules <- names(nHighModule)[nHighModule >= nModuleBound]
  
  ## modules.order <- c(globalModules, setdiff(modules.order, globalModules))
  cols.order <- unlist(lapply(modules.order, function(m) {
    names(g)[g == as.integer(gsub("g", "", m))]
  }))
  return(cols.order)
}


plotCorDensity <- function(corRealFile, corRdmFile,
                           pdcCorFile, corMethod = "SCC",
                           outdir = ".",
                           prefix = "sa2.") {
  # * draw correlation and thershold
  corReal <- loadCor(corRealFile)
  corReal <- corReal[!is.na(corReal$cor), ]

  # * plot line and tail area
  corReal.fdr <- data.table::fread(file = pdcCorFile,
    header = TRUE, sep = "\t",
    data.table = FALSE)
  (minpcc <- min(
    with(corReal.fdr, pcc[(pcc >= 0) & (FDR <=0.01) ])))
  p.tmp <- ggplot(data = corReal, aes(x = cor)) +
    geom_density()
  d <- ggplot_build(p.tmp)$data[[1]]
  # * text
  n_pos <- sum(with(corReal.fdr, (pcc >=0) & (FDR <= 0.01)))
  n_neg <- sum(with(corReal.fdr, (pcc < 0) & (FDR <= 0.01)))
  x_text <- minpcc + 0.01
  y_text <- 3
  textPlus <- paste(format(n_pos, big.mark = ",", scientific = FALSE),
    "positively correlated pairs")


  corRdmShuf <- loadCor(corRdmFile)
  corRdmShuf <- corRdmShuf[!is.na(corRdmShuf$cor), ]

  set.seed(2023)
  corRdmShuf <- corRdmShuf[sample(nrow(corRdmShuf), nrow(corReal)),]
  corsum <- rbind(corReal,corRdmShuf)
  corsum$class <- c(
    rep("Real", nrow(corReal)), rep("RdmShuf", nrow(corReal)))

  withr::with_pdf(
    new = file.path(outdir,
      paste0(prefix,"subclass.cor.", corMethod, ".density.pdf")),
    code = {
      p.density <- ggplot() +
        geom_density(data = corRdmShuf, aes(x= cor),
          fill = "lightgray", colour = NA) +
        geom_density(data = corReal, aes(x = cor),
          colour = "black") +
        geom_vline(xintercept = minpcc, linetype = "dashed") +
        geom_area(data = subset(d, x >= minpcc), aes(x = x, y = y),
          fill = "darkred", alpha = 0.4) +
        geom_text(mapping = aes(label = textPlus, x = x_text + 0.2,
          y = y_text, angle = 90), size = 8, family = "serif") + 
        scale_x_continuous(limits = c(-1.05,1.05), expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        xlab(corMethod) +
        ylab("Density") + 
        theme_classic() +
        theme(
          axis.text = element_text(colour = "black", family = "serif", size = 14),
          axis.title = element_text(colour = "black", family = "serif", size =18))
      print(p.density)
    }, width = 7, height = 7)
}

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

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
