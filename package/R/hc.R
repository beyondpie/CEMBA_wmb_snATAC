sort_hclust <- function(...) {
  as.hclust(
    dendsort::dendsort(as.dendrogram(...))
  )
}

plot_hc <- function(hc, cex = 0.8, mar = c(8,8,4,1)) {
  withr::with_par(new = list(
    cex = cex, mar=mar),
    code = {
      plot(as.dendrogram(hc))
    })
}

#' hiearchical clustering for subclass
#' cpm: cell by features
#' cov.cpm: coefficicent of variation for features
#' postcpm: cell by features
#' @export
hc.subclass.default <- function(cpm,
                                cov.cpm,
                                postcpm,
                                low.q.cpm = 0.05,
                                high.q.cpm = 0.95,
                                distance = 'pearson') {
  idx1 <- which(cov.cpm > quantile(cov.cpm, low.q.cpm))
  idx2 <- which(cov.cpm < quantile(cov.cpm, high.q.cpm))
  idx <- sort(intersect(idx1, idx2))
  mat.hc <- postcpm[ , idx]
  dist.hc <- factoextra::get_dist(x = mat.hc, method = distance)
  hc <- hclust(dist.hc, method = "ward.D2")
  return(list(mat = mat.hc, dist = dist.hc, hc = hc))
}

#' @export
getNMFHeatMap <- function(cpm.plot,
                          ha_row = NULL,
                          ha_col = NULL,
                          showRowNames = TRUE,
                          showColNames = FALSE,
                          fontsize = 6,
                          low.val.col = quantile(cpm.plot, 0.01),
                          high.val.col = quantile(cpm.plot, 0.99),
                          use_raster = TRUE,
                          legend_title = "log2 of CPM",
                          legend_labels = c("Low", "High")) {
  col_fun <- circlize::colorRamp2(
    seq(low.val.col, high.val.col, length = 60),
    viridis::viridis(60)
  )
  p <- ComplexHeatmap::Heatmap(
    matrix = cpm.plot,
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = showRowNames,
    row_names_gp = grid::gpar(fontsize = fontsize),
    show_column_names = showColNames,
    top_annotation = ha_col,
    left_annotation = ha_row,
    use_raster = use_raster,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = legend_title,
      at = c(low.val.col, high.val.col),
      labels = legend_labels)
  )
  return(p)
}

#' @export
gethc.classlevel <- function(dist.mat,
                             sc2mc,
                             nn2neuron = 8.0,
                             gaba2glut = 2.0) {
  old.dist.mat <- dist.mat
  mcOfSc <- sc2mc[rownames(dist.mat), 2]
  index.nn <- grep("NN", mcOfSc, ignore.case = TRUE)
  index.glut <- grep("GLUT", mcOfSc, ignore.case = TRUE)
  index.gaba <- grep("GABA", mcOfSc, ignore.case = TRUE)
  dist.mat[index.nn, c(index.glut, index.gaba)] <- nn2neuron
  dist.mat[c(index.glut, index.gaba), index.nn] <- nn2neuron
  dist.mat[index.glut, index.gaba] <- gaba2glut
  dist.mat[index.gaba, index.glut] <- gaba2glut
  dist.mat.class <- dist.mat
  hc <- hclust(as.dist(dist.mat), method = "ward.D2")
  return(list(
    dist.mat = dist.mat.class,
    old.mat.dist = old.dist.mat,
    hc = hc,
    nn2neuron = nn2neuron,
    gaba2glut = gaba2glut
    ))
}
