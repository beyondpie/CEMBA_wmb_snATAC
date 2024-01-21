library(data.table)
library(tidyverse)

# read ours
my_md5sum <- data.table::fread(file = "raw_bam_file.md5sum",
  header = FALSE, data.table = FALSE) |>
  setNames(object = _, nm = c("md5sum", "bam"))
my_md5sum$bam <- gsub("raw_bams/", "", my_md5sum$bam)
rownames(my_md5sum) <- my_md5sum$bam

# read remotely downloaded ones.
remote_md5sum <- data.table::fread(file = "catlas2023_bam_md5sums.txt",
  header = FALSE, data.table = FALSE) |>
  setNames(object = _, nm = c("bam", "md5sum"))
rownames(remote_md5sum) <- remote_md5sum$bam

all(remote_md5sum[my_md5sum$bam, "md5sum"] == my_md5sum$md5sum)


