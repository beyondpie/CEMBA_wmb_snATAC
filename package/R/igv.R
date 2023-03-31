#' @export
setTrackNode <- function(attrkey,
                         name,
                         attrPath,
                         altColor = "0,0,178",
                         auto.scale = "false",
                         fontSize = 15,
                         height = 20,
                         maximum = 10.0,
                         minimum = 0.0) {
  d <- list()
  attr(d, "baseline") <- "0.0"
  attr(d, "drawBaseline") <- "true"
  attr(d, "flipAxis") <- "false"
  attr(d, "maximum") <- as.character(maximum)
  attr(d, "minimum") <- as.character(minimum)
  attr(d, "type") <- "LINEAR"
  r <- list(DataRange = d)
  attr(r, "altColor") <- altColor
  attr(r, "attributeKey") <- attrkey
  attr(r, "autoScale") <- auto.scale
  attr(r, "clazz") <- "org.broad.igv.track.DataSourceTrack"
  attr(r, "color") <- altColor
  attr(r, "colorScale") <-
    "ContinuousColorScale;0.0;10.0;255,255,255;0,0,178"
  attr(r, "fontSize") <- as.character(fontSize)
  attr(r, "height") <- as.character(height)
  attr(r, "id") <- attrPath
  attr(r, "name") <- name
  attr(r, "renderer") <- "BAR_CHART"
  attr(r, "visible") <- "true"
  attr(r, "windowFunction") <- "mean"
  return(list(track = r))
}

# FIXME: data track order is not followed our suggestion, not sure why.
#' @export
addTrackResource <- function(igv,
                             attrkey,
                             path,
                             name,
                             altColor,
                             autoScale,
                             minimum,
                             maximum,
                             fontSize,
                             height,
                             type) {
  # set a resource node
  resourceNode <- list(Resource = list())
  attr(resourceNode$Resource, "path") <- path
  attr(resourceNode$Resource, "type") <- type
  xml2::xml_add_child(
    .x = xml2::xml_child(igv, "Resources"),
    .value = xml2::as_xml_document(resourceNode)
  )
  # set the corresponding track
  trackNode <- setTrackNode(
    attrkey = attrkey,
    attrPath = path,
    name = name,
    altColor = altColor,
    auto.scale = autoScale,
    fontSize = fontSize,
    height = height,
    maximum = maximum,
    minimum = minimum
  )
  xml2::xml_add_child(
    .x = xml2::xml_child(igv, 2),
    .value = xml2::as_xml_document(trackNode)
  )
  return(igv)
}

#' @export
setGeneList <- function(igv, genes, collapse = "\n") {
  node <- xml2::xml_child(igv, "GeneList")
  xml2::xml_text(node) <- paste(
    genes, collapse = collapse)
  return(igv)
}

#' @export
getRefSeqNode <- function(fontSize = 20, height = 15) {
  r <- list()
  attr(r, "attributeKey") <- "Refseq"
  attr(r, "clazz") <- "org.broad.igv.track.FeatureTrack"
  attr(r, "featureNameProperty") <- "Snap25"
  attr(r, "fontSize") <- as.character(fontSize)
  attr(r, "groupByStrand") <- "false"
  attr(r, "height") <- as.character(height)
  attr(r, "id") <- "https://s3.amazonaws.com/igv.org.genomes/mm10/ncbiRefSeq.sorted.txt.gz"
  attr(r, "name") <- "Refseq"
  attr(r, "visible") <- "true"
  refGenomeNode <- list(Track = r)
  attr(refGenomeNode, "height") <- "1051"
  attr(refGenomeNode, "name") <- "DataPanel"
  attr(refGenomeNode, "width") <- "2497"
  return(refGenomeNode)
}
