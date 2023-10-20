library(tidyverse)
projdir <- here::here()
rpack <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rpack,
  Sa2PeakCalling, cembav2env)
globalPeaks <- Sa2PeakCalling$readBed4File(
  file.path(projdir, "20.nmf", "out",
    "peaks.from.global.modules.bed")
)

pdc <- data.table::fread(
  file = cembav2env$sa2.all.pdc.info, header = TRUE, sep = "\t",
  data.table = FALSE)

global.pdc <- pdc[pdc$cre2 %in% globalPeaks$name, ] |>
  group_by(cre1, cre2) |>
  slice_max(coaccess, n = 1)

global.pdc.srt <- global.pdc |> arrange(desc(coaccess))

proximal <- vapply(global.pdc.srt$cre1,
  str_split_1, FUN.VALUE = rep("1",3), pattern = ":|-") |>
  t() |>
  as.data.frame()

distal <- vapply(global.pdc.srt$cre2,
  str_split_1, FUN.VALUE = rep("1",3), pattern = ":|-") |>
  t() |>
  as.data.frame()

global.pdc.srt.bedpe <- cbind(proximal, distal)
global.pdc.srt.bedpe$coaccess <- global.pdc.srt$coaccess
global.pdc.srt.bedpe$gene <- global.pdc.srt$gene

write.table(global.pdc.srt.bedpe,
  file = file.path(projdir, "04.cCREgene", "sa2.cicero",
    "out/sa2pdcsum", "mba.whole.sa2subclass.globalbynmf.pdc.bedpe"),
  quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
