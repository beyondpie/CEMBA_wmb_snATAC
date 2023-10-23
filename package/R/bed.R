#' @export
loadGRFromBed <- function(bedFile,
                           header = FALSE,
                           sep = "\t",
                           colnms = c("chr", "start", "end", "name")){
  beds <- data.table::fread(file = bedFile,
    sep = sep, header = header)
  if(!header) {
    colnms <- colnms
    colnames(beds) <- colnms
  } else {
    colnms <- colnames(beds)
  }
  gr <- GenomicRanges::GRanges(
    seqnames = beds[[1]],
    ranges = IRanges::IRanges(
      start = beds[[2]],
      end = beds[[3]])
  )
  for(j in seq(from = 4, to = length(colnms))){
    S4Vectors::mcols(gr)[, colnms[j]] <- beds[[j]]
  }
  return(gr)
}

#' @export
strGRSeq <- function(gr) {
  chrs.tmp <- GenomicRanges::seqnames(gr)
  chrs.tmp <- as.factor(chrs.tmp)
  chrs <- levels(chrs.tmp)[chrs.tmp]
  starts <- GenomicRanges::start(gr)
  ends <- GenomicRanges::end(gr)
  r <- paste(chrs, paste(starts, ends, sep = "-"), sep = ":")
  return(r)
}

#' Modification of GenomicRanges::findOverlaps by adding
#' overlap region, width.
#' @export
overlapGR1 <- function(queryGR,
                       subjectGR){
  hits <- GenomicRanges::findOverlaps(
    query = queryGR, subject = subjectGR)
  queryIndex <- S4Vectors::queryHits(hits)
  subjectIndex <- S4Vectors::subjectHits(hits)
  ## get intersect regions
  s <- GenomicRanges::pintersect(
    x = queryGR[queryIndex],
    y = subjectGR[subjectIndex])
  s$ovlpWidth <- GenomicRanges::width(s)
  ## add subject and query seqs
  s$queryGR <- queryGR[queryIndex]
  s$subjectGR <- subjectGR[subjectIndex]
  ## restore hits info
  s$hit <- NULL # remove original TRUE label
  s$ovlp <- hits
  return(s)
}

#' @export
outGR2Bed1 <- function(gr, zeroBased = TRUE) {
  outdf <- as.data.frame(outGR2)
  for( i in colnames(outdf)) {
    if(is.factor(outdf[,i])) {
      message("Transform factor to characters: ", i)
      outdf[, i] <-  levels(outdf[, i])[outdf[,i]]
    }
  }
  outdf$strand <- gsub("\\*", ".", outdf$strand)
}

#' @export
getGeneGR <- function(gtf,
                       mcolnms = c("gene_id", "type", "gene_type", "gene_name")) {
  types <- S4Vectors::mcols(gtf)$type
  types <- levels(types)[types]
  index.gene <- types %in% 'gene'
  geneGTF <- gtf[index.gene, mcolnms]
  return(geneGTF)
}

#' @export
aroundGeneTSS <- function(gtf,
                      up = 1000,
                      down = 1000,
                      mcolnms = c("gene_id", "type", "gene_type", "gene_name"),
                      ...) {
  types <- S4Vectors::mcols(gtf)$type
  types <- levels(types)[types]
  index.gene <- types %in% 'gene'
  geneGTF <- gtf[index.gene, mcolnms]
  r <- GenomicRanges::promoters(
    geneGTF, upstream = up, downstream = down, ...)
  return(r)
}

#' @export
getGeneTSS <- function(gtf,
                       mcolnms = c("gene_id", "type", "gene_type", "gene_name"),
                       ...)  {
  r <- aroundGeneTSS(gtf = gtf, up = 0, down = 1, mcolnms = mcolnms, ...)
  return(r)
}

#' @export
resizeGeneTSS <- function(gtf,
                      window.size = 1000,
                      mcolnms = c("gene_id", "type", "gene_type", "gene_name"),
                      fix = "center",
                      ...) {
  tss <- getGeneTSS(gtf = gtf, mcolnms = mcolnms)
  r <- GenomicRanges::resize(x = tss, width = window.size, fix = fix, ...)
  return(r)
}
