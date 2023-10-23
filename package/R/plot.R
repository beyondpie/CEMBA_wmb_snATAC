library(ggplot2)
library(tidyverse)

# continous color function (from white to red)
colfn <- circlize::colorRamp2(breaks = c(0,1), c("white", "red"))

#' @param group.order: order the group
#' @export
drawReplicateRatio <- function(g2r, group.order = NULL) {
  if (!is.null(group.order)) {
    g2r <- g2r[g2r[, 1] %in% group.order, ]
  }
  r.tbl <- table(g2r[, 1], g2r[, 2])
  r <- as.data.frame(r.tbl / rowSums(r.tbl))
  colnames(r) <- c("group", "sample", "freq")
  if (!is.null(group.order)) {
    r$group <- factor(
      levels(r$group)[r$group], levels = group.order)
  }
  p <- ggplot(r, aes(y = group, x = freq, fill = sample)) +
    geom_bar(position = "fill", stat = "identity")
  return(p)
}

#' @param groupType which column to filter meta
#' @param group keys for filtering meta
#' @param statCol which column to stat
#' @param plotLevelCol which column stat on
#' It can be the same groupType or another groupType
#' @export
cembaBarStat <- function(cemball,
                         groupType = "L2",
                         group = c("NonN_1"),
                         groupOrder = NULL,
                         statCol = "MajorRegion",
                         plotLevelCol = "L3",
                         threshold = 0.0) {
  localM <- cemball[cemball[, groupType] %in% group, ]
  cId2region <- table(localM[, plotLevelCol], localM[, statCol])
  cId2region <- cId2region / rowSums(cId2region)
  cId2region <- reshape2::melt(cId2region)
  colnames(cId2region) <- c("cluster", "region", "ratio")
  if(!is.null(groupOrder)) {
    cId2region$cluster <- factor(cId2region$cluster, levels = groupOrder)
  }
  cId2region <- cId2region[cId2region$ratio > threshold, ]
  return(cId2region)
}

#' @export
cembaBarPlot <- function(cemball,
                         groupType = "L2",
                         group = "NonN_1",
                         groupOrder = NULL,
                         plotLevelCol = "L3",
                         statCol = "MajorRegion",
                         colorCol = NULL,
                         globalTextSize = 20,
                         xAxisAngle = 0,
                         threshold = 0.0) {
  localM <- cemball[cemball[, groupType] %in% group, ]
  cId2region <- table(localM[, plotLevelCol], localM[, statCol])
  cId2region <- cId2region / rowSums(cId2region)
  cId2region <- reshape2::melt(cId2region)
  colnames(cId2region) <- c("cluster", "region", "ratio")
  if(!is.null(groupOrder)){
    cId2region$cluster <- factor(cId2region$cluster, levels = groupOrder)
  }
  cId2region <- cId2region[cId2region$ratio > threshold, ]
  bp <- ggplot(data = cId2region,
    aes(x = cluster, y = ratio, fill = region)) +
    geom_bar(stat = "identity") +
    theme(text = element_text(size = globalTextSize),
      axis.text.x = element_text(angle = xAxisAngle))
  if (!is.null(colorCol)) {
    pal <- cemball[, colorCol]
    names(pal) <- cemball[, statCol]
    uRegion <- unique(names(pal))
    uPal <- pal[uRegion]
    bp <- bp + scale_fill_manual(values = uPal)
  }
  return(bp)
}

#' @import ggplot2
#' @export
vlnPlot <- function(mat, markers,
                    pals = NULL,
                    ylab = "log_CPM",
                    trim = TRUE) {
  sMarkers <- markers[markers %in% colnames(mat)]
  vlnmat <- as.matrix(mat[, sMarkers])
  if(ncol(vlnmat) < length(markers)) {
    message(stringr::str_glue("{ncol(vlnmat)} out of {length(markers)}",
      " markers found "))
  }
  ## to three columns in order to use ggplot2
  vlnmat <- reshape2::melt(vlnmat)
  colnames(vlnmat) <- c("barcode", "gene", ylab)
  vlnmat$barcode <- as.character(vlnmat$barcode)
  vlnmat$gene <- as.character(vlnmat$gene)
  vlnmat$gene <- factor(vlnmat$gene, levels = sMarkers)
  if(is.null(pals)) {
    pals <- grDevices::colorRampPalette(
      pureRUtils::ArchRPalettes$stallion
    )(length(sMarkers))
  }
  p <- ggplot(vlnmat, aes_string(x = "gene", y = ylab, fill = "gene")) +
    geom_violin(trim = trim) +
    scale_fill_manual(values = pals) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20, angle = 90,
        colour = "black"),
      axis.title = element_text(size = 20, colour = "black"),
      legend.position = "none")
  return(p)
}

#' @avg.exp.scaled: gene by ids
#' @pct.exp: gene by ids
#' @import ggplot2
#' @import Seurat
#' @export
mydotplot <- function(avg.exp.scaled,
                      pct.exp,
                      features,
                      ids = NULL,
                      col.min = -2.5,
                      col.max = 2.5,
                      dot.min = 0,
                      dot.scale = 6,
                      colors = c("lightgrey", "blue"),
                      scale.by = 'radius',
                      scale.min = NA,
                      scale.max = NA) {
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("scale.by should be radius or size.")
  )
  avgExpScaled <- Seurat::MinMax(
    data = avg.exp.scaled,
    min = col.min,
    max = col.max)
  avgExpScaled <- if(!is.null(ids)) {
    avgExpScaled[features, ids]
    pctExp <- pct.exp[features, ids]
  } else {
    avgExpScaled[features, ]
    pctExp <- pct.exp[features, ]
  }
  pctExp <- as.vector(pctExp)
  pctExp[pctExp < dot.min] <- NA
  pctExp <- pctExp * 100
  if(!is.na(scale.min)) {
    pctExp[pctExp < scale.min] <- scale.min
  }
  if(!is.na(scale.max)) {
    pctExp[pctExp > scale.max] <- scale.max
  }
  
  ids <- rep(colnames(avgExpScaled), each = nrow(avgExpScaled))
  genes <- rep(rownames(avgExpScaled), times = ncol(avgExpScaled))

  ## colorPool <- grDevices::colorRampPalette(colors = colors)(20)
  ## avgExpBreak <- as.numeric(cut(avgExpScaled, breaks = 20))
  ## colorsAvgExp <- colorPool[avgExpBreak]
  ## by.step <- (col.max-col.min) / 20
  
  dataPlot <- data.frame(
    feature.plot  = factor(x = genes, levels = features),
    id  = ids,
    pct.exp = pctExp,
    avg.exp.scaled = as.vector(avgExpScaled)
  )
  p <- ggplot(data = dataPlot,
    mapping = aes_string(
      x = "feature.plot", y = "id")) +
    geom_point(
      mapping = aes_string(
        size = "pct.exp", color = "avg.exp.scaled")) +
    scale.func(range = c(0, dot.scale),
      limits = c(scale.min, scale.max)) +
    scale_color_gradient(
      low = colors[1], high = colors[2],
      ## breaks = seq(from = col.min, to = col.max, by = by.step),
      limits = c(col.min, col.max)) + 
    theme(axis.title.x = element_blank(),
      axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(x = 'Features', y = 'Annotations') +
    cowplot::theme_cowplot()
  return(p)
}

#' @export
focus.ggplot <- function(p, labels, colors,
                         focus = "dopa",
                         color.bg = "gray") {
  index <- grepl(focus, labels, ignore.case = TRUE)
  if (sum(index) < 1) {
    stop("Cannot find focus: ", focus, " in the labels.")
  }
  colors[!index] <- color.bg
  pp <- p + ggplot2::scale_color_manual(
    labels = labels, values = colors
  )
  return(pp)
}

#' Depend on pureRutils::plotFeatureSingle
#' @export
focus.plotFeatureSingle <- function(featureValue,
                                    focus = "dopa",
                                    color.focus = NULL,
                                    color.bg = "gray",
                                    pal = NULL,
                                    label.not.focus = "others",
                                    keep.label = FALSE,
                                    ignore.case = TRUE,
                                    showLegend = TRUE,
                                    title = focus,
                                    colorOrder = NULL,
                                    ...) {
  index <- grepl(focus, featureValue, ignore.case = TRUE)
  if (sum(index) < 1) {
    stop("Cannot find focus: ", focus, " in the featureValue.")
  }
  labels.match.focus <- unique(featureValue[index])
  newFeatureValue <- featureValue
  if (!keep.label) {
    newFeatureValue[!index] <- label.not.focus
    colorOrder <- c(labels.match.focus, label.not.focus)
    if (is.null(color.focus)) {
      if (is.null(pal)) {
        stop("color.focus and pal are NULL at the same time.")
      }
      color.focus <- pal[labels.match.focus]
    } else {
      if (length(color.focus) < length(labels.match.focus)) {
        warning("Match more labels than the",
          "number of the given color.focus.")
        color.focus <- rep(color.focus[1], length(labels.match.focus))
      }
    }
    pal <- c(color.focus, color.bg)
    names(pal) <- colorOrder
  } else {
    ## set pal
    uLabels <- unique(featureValue)
    if (is.null(pal)) {
      pal <- rep(color.bg, uLabels)
      names(pal) <- uLabels
      if (is.null(color.focus)) {
        stop("No pal and color.focus at the same time.")
      }
      pal[!(names(pal) %in% labels.match.focus)] <- color.bg
    } else {
      if (!all(uLabels %in% names(pal))) {
        stop("Not all the labels have the corresponding colors in pal.")
      }
    }
    if (!is.null(color.focus)) {
      pal[labels.match.focus] <- color.focus
    }
  }
  p <- pureRUtils::plotFeatureSingle(
    featureValue = newFeatureValue,
    showLegend = showLegend,
    pal = pal,
    colorOrder = colorOrder,
    ...
  ) + ggplot2::ggtitle(title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5))
  return(p)
}


#' @export
# https://github.com/yanlinlin82/ggvenn/issues/13
gen_circle <- function(group, x_offset = 0, y_offset = 0, radius = 1,
                       radius_b = radius, theta_offset = 0, length.out = 100) {
  tibble(group = group,
         theta = seq(0, 2 * pi, length.out = length.out)) %>%
    mutate(x_raw = radius * cos(theta),
           y_raw = radius_b * sin(theta),
           x = x_offset + x_raw * cos(theta_offset) - y_raw * sin(theta_offset),
           y = y_offset + x_raw * sin(theta_offset) + y_raw * cos(theta_offset))
}

#' @export
prepareVennData.d3 <- function(venn.data) {
  label.text <- names(venn.data)
  size.each <- sapply(venn.data, length)
  size.intersect <- c(
    length(intersect(venn.data[[1]], venn.data[[2]])),
    length(intersect(venn.data[[1]], venn.data[[3]])),
    length(intersect(venn.data[[2]], venn.data[[3]])),
    length(intersect(intersect(venn.data[[1]], venn.data[[2]]),
      venn.data[[3]]))
    )
  size.intersect[1:3] <- size.intersect[1:3]- size.intersect[4]
  size.each2 <- size.each
  size.each2[1] <- length(setdiff(setdiff(venn.data[[1]], venn.data[[2]]),
    venn.data[[3]]))
  size.each2[2] <- length(setdiff(setdiff(venn.data[[2]], venn.data[[1]]),
    venn.data[[3]]))
  size.each2[3] <- length(setdiff(setdiff(venn.data[[3]], venn.data[[1]]),
    venn.data[[2]]))
  text.venn <- c(size.each2, size.intersect)
  names(text.venn) <- c("x1", "x2", "x3", "x12", "x13", "x23", "x123")
  return(list(label = label.text, text = text.venn))
}

#' @export
drawVenn3 <- function(text.to.use,
                      label.to.use,
                      text.size = 15,
                      label.size = 12,
                      radius = 2.0, scale = 0.8,
                      l1 = radius * 0.5, l2 = radius * 0.7,
                      l3 = radius * 0.1,
                      fontface = "bold",
                      edge.colors = c("blue", "darkgray", "darkgreen"),
                      fill.colors = c("red", "yellow", "green")
                      ) {
  x1 <- 0.0
  y1 <- 0.0
  x2 <- x1 + radius * scale
  y2 <- y1
  x3 <- x1 + scale * radius / 2
  y3 <- y1 + scale * radius * sqrt(3)/2

  x11 <- x1 - cos(pi/4) * l1
  y11 <- y1 - cos(pi/4) * l1
  x22 <- x2 + cos(pi/4) * l1
  y22 <- y2 - cos(pi/4) * l1
  x33 <- x3
  y33 <- y3 + l1

  x12 <- (x1 + x2) / 2
  y12 <- y1 - l2
  x13 <- x1
  y13 <- y1 + l2
  x23 <- x2
  y23 <- y2 + l2
  x123 <- x3
  y123 <- y3 - l2

  text.pos <- data.frame(
    x = c(x11, x22, x33, x12, x13, x23, x123),
    y = c(y11, y22, y33, y12, y13, y23, y123),
    text = text.to.use
  )
  label.pos <- data.frame(
    x = c(x1, x2, x3),
    y = c(y1-radius-l3, y2-radius-l3, y3 + radius + l3),
    name = label.to.use
  )
  g <- ggplot() +
    geom_polygon(aes(x = x, y = y),
      data = gen_circle(group = 0L, x_offset = x1, y_offset = y1, radius = radius),
      color = edge.colors[1], fill = fill.colors[1], alpha = .5) +
    geom_polygon(aes(x = x, y = y),
      data = gen_circle(group = 0L, x_offset = x2, y_offset = y2, radius = radius),
      color = edge.colors[2], fill = fill.colors[2], alpha = .5) +
    geom_polygon(aes(x = x, y = y),
      data = gen_circle(group = 0L, x_offset = x3, y_offset = y3, radius = radius),
      color = edge.colors[3] , fill = fill.colors[3], alpha = .5) +
    geom_text(data = text.pos, aes(x = x, y = y, label = text),
      color = "black", size = text.size, fontface = fontface, family = "serif") +
    geom_text(data = label.pos, aes(x = x, y = y, label = name),
      color = "black", size = label.size, family = "serif", fontface = fontface) + 
    coord_fixed(ratio = 1) + 
    theme_void()
  return(g)
}

#' @export
get.dissect.desc <- function(atacMeta) {
  cell.slice.h <- atacMeta$Slice
  cell.slice.v <- gsub("\\d+", "", atacMeta$dissect)
  slice.h <- sort(unique(cell.slice.h))
  slice.v <- sort(unique(cell.slice.v))
  n.slice <- table(cell.slice.v, cell.slice.h)
  size.slice <- as.data.frame(n.slice)
  size.slice$cell.slice.h <- levels(
    size.slice$cell.slice.h)[size.slice$cell.slice.h]
  size.slice$cell.slice.v <- levels(
    size.slice$cell.slice.v)[size.slice$cell.slice.v]
  rownames(size.slice) <- paste0(
    size.slice$cell.slice.h, size.slice$cell.slice.v
  )
  mr.hv <- rep("empty", nrow(size.slice))
  v1.hv <- rep("FALSE", nrow(size.slice))
  for( i in seq_len(nrow(size.slice))) {
    dissect <- rownames(size.slice)[i]
    index <- atacMeta$dissect == dissect
    if (sum(index) > 0) {
      r <- table(atacMeta$MajorRegion[index])
      str.r <- paste(names(r), r, sep = "_", collapse = "-")
      mr.hv[i] <- names(r)[which.max(r)]

      r2 <- table(atacMeta$CEMBAv1[index])
      str.r2 <- paste(names(r2), r2, sep = "_", collapse = "-")
      v1.hv[i] <- names(r2)[which.max(r2)]
      message(paste(dissect, str.r, str.r2, sep = ":"))
    }
  }
  desc <- data.frame(
    sliceH = size.slice$cell.slice.h,
    sliceV = size.slice$cell.slice.v,
    mr = mr.hv,
    v1 = v1.hv,
    size = size.slice$Freq
  )
  return(desc)
}

