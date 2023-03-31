#' @examples
#' File format:
#' pdc,cor
#' 0610009B22Rik@peak217336,0.0602667683204699
#' 0610009B22Rik@peak217343,-0.0929365338771179
#' @export
loadCor <- function(f, split = "@") {
  r <- data.table::fread(file = f, header = TRUE, sep = ",",
    data.table = FALSE)
  gene.peak  <- vapply(r$pdc, function(i) {
    t <- strsplit(i, split = split, fixed = TRUE)
    t[[1]]
  }, c("a", "peak1"))
  r$gene <- gene.peak[1,]
  r$peak <- gene.peak[2,]
  return(r)
}
