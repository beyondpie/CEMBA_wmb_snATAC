library(purrr)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(stringr)
projdir <- here::here()
resourcedir <- "../resource"
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  Sa2Integration)

# * simply extract allen's cluster marker genes
annotAllen <- Sa2Integration$loadAllenAnnot()

# * functions
filterMarkers <- function(markers, sep = ",") {
  # filter space
  # filter gene symbols linked with ENSMUSG
  r <- map(markers, str_split_1, pattern = sep) |>
    unlist() |> sort() |>
    x => x[!grepl("ENSMUSG", x)] |>
    x => x[nzchar(x)] |>
    unique()
  message("found ", length(r), " markers.")
  message("unique ", length(unique(r)), " markers.")
  return(r)
}
getMarkers <- function(annotAllen,
                       groupBy = "",
                       groups = c(),
                       markerCol = "cluster.markers",
                       outfnm = "") {
  message("marker column: ", markerCol)
  markers <- if (nchar(groupBy) > 1) {
    message("select markers from column: ", groupBy)
    message("with values: ", paste(groups, collapse = " "))
    gs <- annotAllen[[groupBy]]
    annotAllen[gs %in% groups, markerCol]
  } else {
    annotAllen[[markerCol]]
  }
  r <- filterMarkers(markers)
  if (nchar(outfnm) > 1) {
    message("save results to: ", outfnm)
    write.table(r, outfnm, append = FALSE, quote = FALSE,
      col.names = FALSE, row.names = FALSE)
  }
  return(r)
}

# * main
allMarkers <- getMarkers(
  annotAllen, markerCol = "cluster.markers",
  outfnm = file.path(resourcedir, "AIT21_cluster_markers.txt")
)
allMerfishMarkers <- getMarkers(
  annotAllen, markerCol = "merfish.markers",
  outfnm = file.path(resourcedir, "AIT21_merfish_markers.txt")
)

# fix k8 markers
rawk8markers <- read.table(
  Sa2Integration$AllenRaw8KMarkerFile, header = FALSE)$V1
k8markers <- filterMarkers(rawk8markers)
write.table(k8markers,
  file.path(resourcedir, "AIT21_k8_markers.txt"),
  append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)

# * test markers
geneList <- Sa2Integration$getMarkersList()
