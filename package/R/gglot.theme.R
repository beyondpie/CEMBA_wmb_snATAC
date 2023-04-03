library(ggplot2)

theme.no.axis <- theme(
  axis.text = ggplot2::element_blank(),
  axis.title = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(),
  panel.border  = ggplot2::element_blank(),
  panel.background = element_blank(),
  plot.margin = unit(c(0,0,0,0), "cm"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 9)
)

theme_my_minimal <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_text(colour = "black"),
    )

#'@export
setGlobalTheme <- function(newTheme){
  oldTheme <- ggplot2::theme_get()
  ggplot2::theme_set(newTheme)
  return(oldTheme)
}
