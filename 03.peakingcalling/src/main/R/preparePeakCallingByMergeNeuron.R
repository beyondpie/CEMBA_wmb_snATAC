library(tidyverse)
library(rlang)
library(R6)

Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
projdir <- here::here()
rdir <- file.path(projdir, "package/R")
import::from(.from = "cembav2env.R", .directory = rdir,
  cembav2env, cluSumBySa2)
import::from(.from = "peak.R", .directory = rdir,
  mapL4pc2L4)

# load tf light summary
fromdir <- file.path(projdir, "19.snap2_integration", "out",
  "tf_lightsum")

# best parameters for transfer label
## kAnchors <- c(30, 50, 70)
kAnchors <- c(50)
method <- "cca"
allenTech <- "allen"
group <- "neuron"
feaName <- "k8"
# use for merging smaller clusters
ncheck <- 200

tfLightSumList <- map(kAnchors, \(k) {
  fnm <- file.path(fromdir,
    str_glue("{group}_atac.{allenTech}_{method}_{feaName}_{k}.lightsum.rds"))
  r <- readRDS(fnm)
  rownames(r) <- r$L4
  colnms <- colnames(r)
  index <- !grepl("L4", colnms)
  colnms[index] <- paste(paste0("k", k), colnms[index], sep = "_")
  colnames(r) <- colnms
  r
})
names(tfLightSumList) <- paste("k", kAnchors, sep = "")

# merge cl_tf, cl_tfscore from the list above
tfLightSum <- Reduce(function(x, y) {
  merge(x, y, by = "L4", all = TRUE)
}, tfLightSumList)
rownames(tfLightSum) <- tfLightSum$L4

# k=50 performs better, same cl and subclass as 70
# k=50 also shows more consistent on subclass with k=30

tfLightSum2 <- tfLightSum[,
  c("L4", "k50_cl_tf", "k50_subclass_id", "L4size")]

#' assume L4s have the same L3.
mergeL4nms <- function(L4s, prefix = "p", link = "_") {
  L3 <- unique(cluSumBySa2$getL3sa2(L4s))
  if (length(L3) > 1) {
    stop("this fn only works within one L3 cluster.")
  }
  L4ids <- gsub("^\\d+-\\d+-\\d+-", "", L4s) |>
    x => paste(x, collapse = link)
  str_glue(
    "{L3}-{prefix}{L4ids}")
}

L4Group4PeakCalling <- R6::R6Class(
  classname = "L4Group4PeakCalling",
  public = list(
    L3 = NULL,
    rawL4s = list(),
    prefix = "p",
    pL4s = list(),
    # cl, subclass, size info
    rawL4Meta = data.frame(),
    initialize = function(L4meta) {
      L4s <- L4meta$L4
      # check they are from the same L3.
      L3s <- cluSumBySa2$getL3sa2(L4s)
      if (length(unique(L3s)) > 1) {
        stop("L4s are not from the same L3.")
      }
      self$L3 <- unique(L3s)
      self$rawL4s <- L4s
      self$rawL4Meta <- L4meta
    },
    mergeL4 = function(groupBy = "k50_cl_tf", newGroupName = "cl") {
      if (length(self$rawL4s) == 1) {
        ## message("only one group in L4, so use it directly.")
        self$pL4s <- data.frame(
          pL4 = self$rawL4s,
          pL4size = self$rawL4Meta$L4size,
          cl = self$rawL4Meta[[groupBy]])
      } else {
        r <- self$rawL4Meta |>
          dplyr::group_by(!!!rlang::syms(groupBy)) |>
          summarize(pL4 = mergeL4nms(L4), pL4size = sum(L4size))
        self$pL4s <- r[, c("pL4", "pL4size", groupBy)] |>
          setNames(object = _, c("pL4", "pL4size", newGroupName))
      }
      invisible(self)
    }
))


smallLightSum2 <- subset(tfLightSum2, L4size < ncheck)
smallLightSum2$L3 <- cluSumBySa2$getL3sa2(smallLightSum2$L4)

mergeLightSumList <- smallLightSum2 |>
  group_split(L3) |>
  map(.x = _, .f = L4Group4PeakCalling$new) |>
  x => map(.x = x, .f = \(y) {
    y$mergeL4(groupBy = "k50_cl_tf", newGroupName = "cl")
  }) |>
  x => setNames(object = x, nm = map_chr(.x = x, \(y) {y$L3}))

pL4s <- map(mergeLightSumList, .f = \(x) x$pL4s) |>
  do.call(what = rbind, args = _) |>
  setNames(object = _, c("L4", "L4size", "cl"))

leftLightSum2 <- subset(tfLightSum2, L4size >= ncheck) |>
  x => x[, c("L4", "L4size", "k50_cl_tf")] |>
  setNames(object = _,  c("L4", "L4size", "cl"))

# *  merge all the L4pc together
L4pc2meta <- rbind(leftLightSum2, pL4s)
rownames(L4pc2meta) <- L4pc2meta$L4

# * attach to barcodes

L4pc2L4 <- mapL4pc2L4(L4pc = L4pc2meta$L4)
rownames(L4pc2L4) <- L4pc2L4$L4

barcode2L4 <- cluSumBySa2$loadbarcode2L4()
b4pc <- barcode2L4[barcode2L4$L4 %in% rownames(L4pc2L4), ]
b4pc['L4pc'] <- L4pc2L4[b4pc$L4, 'pL4']
b4pc['psize'] <- L4pc2meta[b4pc$L4pc, "L4size"]

# add sample rep info to each barcode
sa2Meta <- readRDS(cembav2env$sa2metaFile)
rownames(sa2Meta) <- sa2Meta$barcode2
b4pc['biorep'] <- sa2Meta[b4pc$barcode, "biorep"]
outfnm <- file.path(projdir, "18.snap2_peakcalling",
  "src/main/resource",
  str_glue("{group}_barcode2cluster_bedtag-all.csv"))

write.table(x = b4pc[, c("barcode", "L4pc")],
  file = outfnm, quote= FALSE, row.names = FALSE,
  sep = ",", col.names = FALSE)

# * prepare directly the L4pc_biorep
b4pc$crep <- paste(b4pc$L4pc, b4pc$biorep, sep = ".")
statcrep <- table(b4pc["crep"])
ignorecrep <- names(statcrep)[statcrep < ncheck]
b2crep <- b4pc[!(b4pc$crep %in% ignorecrep), c("barcode", "crep")]
write.table(x = b2crep,
  file = file.path(projdir, "18.snap2_peakcalling",
    "src/main/resource",
    str_glue("{group}_barcode2cluster_bedtag-biorep.csv")),
  quote = FALSE, row.names = FALSE, sep = ",",
  col.names = FALSE)

# * prepare L4pc_pseudorep
getpseudorep <- function(L4pc,
                         repnm = c("early", "later"),
                         prob = c(0.5, 0.5),
                         suffix = 'pseudo',
                         sep = ".") {
  pseudoreps <- sample(repnm, size = length(L4pc),
    replace = TRUE, prob = prob)
  paste(L4pc, pseudoreps, suffix, sep = sep)
}
b2cpseurep <- b4pc |> dplyr::group_by(L4pc) |>
  dplyr::reframe(barcode = barcode, pseudorep = getpseudorep(L4pc)) |>
  as.data.frame()
rownames(b2cpseurep) <- b2cpseurep$barcode

## check result
all(b2cpseurep[b4pc$barcode, "L4pc"] == b4pc$L4pc )
tmp <- map_chr(b2cpseurep$pseudorep, \(x)str_split_1(x, '\\.')[1])
all(tmp ==  b2cpseurep$L4pc)
a <- table(b4pc$L4pc)
cl <- names(a)[100]
(t <- table(b2cpseurep$pseudorep[b2cpseurep$L4pc %in% cl]))
sum(t) == a[100]
## save result
statprep <- table(b2cpseurep$pseudorep)
ignoreprep <- names(statprep)[statprep < ncheck]
b2prep <- b2cpseurep[!(b2cpseurep$pseudorep %in% ignoreprep),
  c("barcode", "pseudorep")]

write.table(x = b2prep,
  file = file.path(projdir, "18.snap2_peakcalling",
    "src/main/resource",
    str_glue("{group}_barcode2cluster_bedtag-pseudorep.csv")),
  quote = FALSE, row.names = FALSE, sep = ",",
  col.names = FALSE)


# * add size, and size of replicates
sizeOfL4pc <- L4pc2meta$L4size
L4pc2rep <- with(b4pc, table(L4pc, biorep))
sizeOfEarlyOfL4pc <- L4pc2rep[L4pc2meta$L4, "early"]
sizeOfLaterOfL4pc <- L4pc2rep[L4pc2meta$L4, "later"]

L4pc2sizes <- data.frame(
  pL4 = L4pc2meta$L4,
  size = sizeOfL4pc,
  sizeEarly = sizeOfEarlyOfL4pc,
  sizeLater = sizeOfLaterOfL4pc
)
write.table(
  x = L4pc2sizes,
  file = file.path(projdir, "18.snap2_peakcalling",
    "src/main/resource",
    str_glue("{group}_L4pc2sizes_{method}-k{kAnchors}-cl_v1.csv")),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = ","
)

# * check some clusters.
l4 <- "10-1-2-1"
l4 <- "10-3-2-1"
samples<- b4pc$barcode[b4pc$L4pc %in% l4] |>
  map_chr(.x = _, .f = \(y) str_split_1(y, '\\.')[1])
table(samples)
