projroot <- here::here()
rdir <- file.path(projroot, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env, Sa2Integration, Sa2PeakCalling)
import::from(.from = "utils.R", .directory = rdir,
  fastread.csv)


allenMeta <- Sa2Integration$loadAllenAnnot()

subclassMap <- unique(data.frame(
  subclass_id = allenMeta$subclass_id,
  subclass_id_label = allenMeta$subclass_id_label,
  subclass_label = allenMeta$subclass_label
))

subclass_label_bw <- data.table::fread(
  file = "../resource/subclass_nm_in_macs2_bigwig.txt",
  header = FALSE, data.table = FALSE)$V1

subclassMap$subclass_id_label_bw <- subclassMap$subclass_id_label |>
  gsub(" ", "_", x = _) |>
  gsub("/", ".", x = _)

subclassMap$subclass_label_bw <- subclassMap$subclass_label |>
  gsub(" ", "_", x = _) |>
  gsub("/", ".", x = _)

subclassMap$subclass_label_peak <- subclassMap$subclass_label |>
  Sa2PeakCalling$renameAllenSubclass()

# check
all(subclass_label_bw %in% subclassMap$subclass_id_label_bw)
a <- fastread.csv(fnm = file.path(projroot, "18.snap2_peakcalling",
  "out/scfilter", "count_peakBysubclass.csv"))
all(colnames(a) %in% subclassMap$subclass_label_peak)

# save
saveRDS(subclassMap,
  file.path(projroot, "meta", "sa2.subclass.names.map.rds"))
write.table(subclassMap,
  file.path(projroot, "meta", "sa2.subclass.names.map.csv"),
  sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
