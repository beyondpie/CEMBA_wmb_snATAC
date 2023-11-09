source("paper.R")
scnmMap <- cembav2env$load.sa2.subclass.names()
rownames(scnmMap) <- scnmMap$subclass_id_label
allenSubclass2Region <- Sa2Integration$getAllenSubclass2Region()
rownames(allenSubclass2Region) <- allenSubclass2Region$subclass_id_label

mrMeta <- list(
  Telencephalon = c("Isocortex", "HPF", "OLF", "AMY", "STR", "PAL"),
  Diencephalon = c("TH", "HY"),
  Midbrain = "MB",
  Hindbrain = c("Pons", "MY"),
  Cerebellum = "CB"
)
mr2r <- data.frame(
  mr = c(
    c("Isocortex", "HPF", "OLF", "AMY", "STR", "PAL", "CNU"),
    c("TH", "HY"),
    c("MB", "MB-PONS"),
    c("Pons", "MY", "HB"),
    "CB"),
  r = c(rep("Telencephalon", 7),
    rep("Diencephalon", 2),
    rep("Midbrain", 2),
    rep("Hindbrain", 3),
    "Cerebellum")
)
rownames(mr2r) <- mr2r$mr

r2Color <- data.frame(
  r = c("Telencephalon", "Diencephalon", "Midbrain",
    "Hindbrain", "Cerebellum", "Non-Telencephalon"),
  color = c("#00688B", "#F15F30", "#74E44B",
    "#788FC8", "#DEB34C", "#999999")
)

# use this to update mr2r
allenmrs <- allenSubclass2Region[, "MajorRegionRateTop3"] |>
  lapply(X = _, Sa2Integration$extractAllenRegionWithScore) |>
  do.call(rbind, args = _ )|>
  x => x[!is.na(x[,2]), ]|>
  x => x$region |>
  unique()

# map subclass to region based on their top region information
sc2r <- allenSubclass2Region[ , "MajorRegionRateTop3"] |>
  lapply(X = _, \(x) {
    r2score <- Sa2Integration$extractAllenRegionWithScore(x) |>
      mutate(r = mr2r[region, "r"]) |>
      group_by(r) |>
      summarise(r_s = sum(score)) |>
      slice_max(r_s, n = 1)
  }) |> do.call(rbind, args = _) |>
  as.data.frame()

rownames(scnmMap) <- scnmMap$subclass_id_label
rownames(sc2r) <- scnmMap[allenSubclass2Region$subclass_id_label,
  "subclass_label_peak"]

write.table(sc2r,
  file = file.path(projdir, "meta", "sa2.subclass2region2score.csv"),
  sep = ",",
  col.names = TRUE, row.names = TRUE, quote = FALSE)

