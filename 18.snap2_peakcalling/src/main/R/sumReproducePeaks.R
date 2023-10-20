## library(R.utils)
## library(purrr)
args <- commandArgs(trailingOnly = TRUE)
rpdpeakDir <- args[1]
outdir <- args[2]

# debug
## rpdpeakDir <-  "/oasis/tscc/scratch/szu/projects/CEMBA2/18.snap2_peakcalling/out/tscc/rpdpeak"
## outdir <- "."

recordSummitFiles <- file.path(outdir, "mba.whole.naiveSummitList.list")
## cl2npeakFile <- file.path(outdir, "mba.whole.L4.npeak4anno.txt")

summitFiles <- list.files(path = rpdpeakDir, pattern = ".*naiveSummitList.bed",
  full.names = TRUE, include.dirs = TRUE, no.. = TRUE, ignore.case = FALSE)

cls <- basename(summitFiles) |> gsub(".naiveSummitList.bed", "", x = _)
## npeaks <- map_int(.x = summitFiles, .f = countLines)

r1 <- data.frame(
  cl = cls,
  path = summitFiles
)

write.table(r1, file = recordSummitFiles, quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = FALSE, append = FALSE)

## r2 <- data.frame(
##   cl = cls,
##   npeak = npeaks
## )
## write.table(r2, file = cl2npeakFile, quote = FALSE, sep = "\t",
##   row.names = FALSE, col.names = FALSE, append = FALSE)

message("done.")

