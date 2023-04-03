# ref: https://rpubs.com/TX-YXL/662586
#' @export
dendro_data_k <- function(hc, k) {
  hcdata <- ggdendro::dendro_data(hc, type = "rectangle")
  seg <- hcdata$segments
  labclust <- cutree(hc, k)[hc$order]
  segclust <- rep(0L, nrow(seg))
  heights <- sort(hc$height, decreasing = TRUE)
  height <- mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)

  for (i in 1:k) {
    xi <- hcdata$labels$x[labclust == i]
    idx1 <- seg$x >= min(xi) & seg$x <= max(xi)
    idx2 <- seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3 <- seg$yend < height
    idx <- idx1 & idx2 & idx3
    segclust[idx] <- i
  }

  idx <- which(segclust == 0L)
  segclust[idx] <- segclust[idx + 1L]
  hcdata$segments$clust <- segclust
  hcdata$segments$line <- as.integer(segclust < 1L)
  hcdata$labels$clust <- labclust
  hcdata
}

#' @export
set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan = FALSE) {
  if (fan) {
    angle <- 360 / nbLabels * 1:nbLabels + 90
    idx <- angle >= 90 & angle <= 270
    angle[idx] <- angle[idx] + 180
    hjust <- rep(0, nbLabels)
    hjust[idx] <- 1
  } else {
    angle <- rep(0, nbLabels)
    hjust <- 0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}

#' @export
plot_ggdendro <- function(hcdata,
                          direction = c("lr", "rl", "tb", "bt"),
                          fan = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size = 3,
                          nudge.label = 0.01,
                          expand.y = 0.1) {
  library(ggplot2)
  library(ggdendro)
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks <- pretty(ggdendro::segment(hcdata)$y, n = 5)
  ymax <- max(ggdendro::segment(hcdata)$y)

  ## branches
  p <- ggplot() +
    geom_segment(data = segment(hcdata),
      aes(x = x,
        y = y,
        ## linetype = factor(line),
        ## colour = factor(label, levels(label))
        xend = xend,
        yend = yend,
      ),
      lineend = "round",
      show.legend = FALSE,
      size = branch.size)
  

  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
        limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip() 
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }

  # labels
  labelParams <- set_labels_params(
    nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle

  p <- p +
    geom_text(data = ggdendro::label(hcdata),
      aes(x = x,
        y = y,
        label = label,
        ## colour = label,
        angle = angle),
      vjust = labelParams$vjust,
      hjust = labelParams$hjust,
      nudge_y = ymax * nudge.label,
      size = label.size,
      show.legend = FALSE,
      color = "black")

  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }

  ylim <- -round(ymax * expand.y, 1)
  p <- p + expand_limits(y = ylim)
  p
}
