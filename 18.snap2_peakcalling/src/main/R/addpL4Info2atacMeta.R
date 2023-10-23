library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env)

atacMeta <- readRDS(cembav2env$sa2metaFile)
rownames(atacMeta) <- atacMeta$barcode2

barcode2pL4_nn_fnm <- file.path(projdir, "18.snap2_peakcalling",
  "src/main/resource",
  "nn_barcode2cluster_bedtag-all.csv")

barcode2pL4_nn <- data.table::fread(
  file = barcode2pL4_nn_fnm, header = FALSE, sep = ",",
  data.table = FALSE)
colnames(barcode2pL4_nn) <- c("barcode", "pL4")
barcode2pL4_nn$gpL4 <- paste0("nn.", barcode2pL4_nn$pL4)

barcode2pL4_neuron_fnm <- file.path(projdir, "18.snap2_peakcalling",
  "src/main/resource",
  "neuron_barcode2cluster_bedtag-all.csv")
barcode2pL4_neuron <- data.table::fread(
  file = barcode2pL4_neuron_fnm, header = FALSE, sep = ",",
  data.table = FALSE)
colnames(barcode2pL4_neuron) <- c("barcode", "pL4")
barcode2pL4_neuron$gpL4 <- paste0("neuron.", barcode2pL4_neuron$pL4)

barcode2pL4 <- rbind(barcode2pL4_nn, barcode2pL4_neuron)
rownames(barcode2pL4) <- barcode2pL4$barcode

atacMeta$pL4 <- barcode2pL4[rownames(atacMeta), "pL4"]
atacMeta$gpL4 <- barcode2pL4[rownames(atacMeta), "gpL4"]
saveRDS(atacMeta,
  file = file.path(projdir, "supple.02.annotation.all.in.one.meta",
    "mba.whole.cell.meta.v9.4.rds"))


# * now it's v9.7
atacMeta <- readRDS(cembav2env$sa2metaFile)

table(atacMeta[grep("IMN", atacMeta$subclass_label_v3), c("L4", "NT_v3")], useNA = "ifany")
table(atacMeta[grep("IMN", atacMeta$subclass_label_v3), c("L4", "mainclass.rough")], useNA = "ifany")
