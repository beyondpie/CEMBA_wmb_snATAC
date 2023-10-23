# * plot entropy vs sparseness
# will detect all the modules in the given directory.

library(ggplot2)
library(optparse)
packdir <- file.path(here::here(), "package/R")
import::from(.from = "utils.R", .directory = packdir,
  checkArgsExistOrStop, prepareOutdir)

op <- list(
  make_option(c("--nmfDir"), type = "character",
    default = "nmf_ppdc/out"),
  make_option(c("--tag"), type = "character",
    default = "ppdc")
)

args <- parse_args(OptionParser(option_list = op))

if(!dir.exists(args$nmfDir)) {
  stop(args$nmfDir, " does not exist.")
}

# * configs
staLiteHeader <- c("cellSparse", "entropy")
staBoxLiteHeader <- c("rank", "Min", "Q1", "Median", "Mean", "Q3",
  "Max", "TopWhisker", "BottomWhisker", "Box1", "Box2", "Box3",
  "UpWhisker", "DnWhisker")


# * function
matchFiles <- function(nmfDir = "nmf_ppdc/out",
                       suffix = "sta.txt$") {
  files <- list.files(path = nmfDir, include.dirs = FALSE,
    no.. = TRUE, recursive = FALSE, pattern = suffix)
  if (length(files) < 1) {
    stop("No matched files under ", nmfDir, " for suffix ", suffix)
  }
  return(files)
}

getNMFPrefix <- function(files, extrastr = ".sta.txt$") {
  v <- vapply(files, gsub,
    pattern = extrastr, replacement = "", FUN.VALUE = "a")
  return(v)
}

getNMFModules <- function(files = NULL,
                          nmfPrefix = NULL,
                          split = ".",
                          index = 3,
                          extrastr = "r") {
  strings <- if(!is.null(files)) {
    files
  } else {
    nmfPrefix
  }
  if(is.null(strings)) {
    stop("No strings for getting nmf modules.")
  }
  mods <- vapply(strings, function(f){
    t <- gsub(extrastr, "",
      unlist(strsplit(f, split = split, fixed = TRUE))[index])
    as.numeric(t)
  }, FUN.VALUE = 1)
  return(mods)
}

# * get nmfPmat.{tag}.r{r}.box.[contributes, sparseness, entropy]
readBoxSta <- function(nmfDir, nmfPrefix, prefix2 = "n0") {
  modules <- getNMFModules(nmfPrefix = nmfPrefix, files = NULL)
  boxstaList <- lapply(seq_along(modules), function(i) {
    boxsta <- read.table(
      file = file.path(nmfDir,
        paste(nmfPrefix[i], prefix2, "box.sta",sep = ".")),
      sep= "\t", header = TRUE, row.names = 1, quote = ""
    )
    return(boxsta)
  })
  chars <- colnames(boxstaList[[1]])
  sumBoxStaList <- lapply(chars, function(nm) {
    r <- do.call(cbind, lapply(boxstaList, function(i) {i[,nm]}))
    rownames(r) <- rownames(boxstaList[[1]])
    colnames(r) <- modules
    return(r)
  } )
  names(sumBoxStaList) <- chars
  return(sumBoxStaList)
}

# * get nmfPmat.{tag}.sta.txt content in data.frame
readStatxt <- function(nmfDir, nmfPrefix, samples = "nmfPmat") {
  modules <- getNMFModules(nmfPrefix = nmfPrefix, files = NULL)
  statxtList <- lapply(seq_along(modules), function(i) {
    s <- data.table::fread(file.path(nmfDir,
      paste0(nmfPrefix[i], ".sta.txt")), header = FALSE,
      data.table = FALSE)$V1
    r <- data.frame(
      samples = rep(samples, 2),
      ranks = rep(modules[i], 2),
      stat = staLiteHeader,
      val = s)
    return(r)
  })
  statxt <- do.call(rbind, statxtList)
  return(statxt)
}

# * get nmfPmat.{tag}.statH.sta.txt  into one data.frame
readStatH <- function(nmfDir, nmfPrefix, prefix2 = "n0",
                      samples = "nmfPmat") {
  modules <- getNMFModules(nmfPrefix = nmfPrefix, files = NULL)
  statHList <- lapply(seq_along(modules), function(i) {
    ## first row has a # char, should be
    ## automatically ignored by read.table
    s <- read.table(file.path(nmfDir,
      paste(nmfPrefix[i], prefix2, "statH", sep = ".")), sep = "\t",
      header = FALSE, quote = "")
    colnames(s) <- c("xgi", "index", "class0",
      "class1", "contributes", "sparseness", "entropy")
    s0 <- data.frame(samples = rep(samples,nrow(s)),
      ranks = rep(modules[i], nrow(s)))
    r <- cbind(s0, s)
    return(r)
  })
  statH <- do.call(rbind, statHList)
  return(statH)
}

plot_entry2sparse <- function(boxStaList,
                              measure = "Median",
                              tag = "ppdc",
                              width = 10,
                              height = 7) {
  entr_stats <- boxStaList$entropy
  sparse_stats <- boxStaList$sparseness
  module_nms <- paste("Module", colnames(entr_stats))
  
  entr_values <- entr_stats[measure, ]
  sparse_values <- sparse_stats[measure, ]
  df_plot <- rbind(entr_values, sparse_values)
  rownames(df_plot) <- c("Entropy", "Sparseness")
  colnames(df_plot) <- module_nms

  mat_plot <- as.matrix(t(df_plot))
  cols <- c("blue", "red")
  withr::with_pdf(
    new = file.path(sum_dir, paste(tag2, "entropy_sparse", "pdf",sep = ".")),
    code = {
      matplot(mat_plot, type = c("b"), pch = 19, col = cols,  xaxt = "n",
        xlab = paste("Number of Modules for", tag),
        ylab = "Measurements",
        main = paste("NMF for", tag))
      axis(side = 1, at = seq_len(nrow(mat_plot)), labels = rownames(mat_plot))
      legend("topright", legend = c("Entropy", "Sparseness"), col = cols, pch = 19)
    }, width = width, height = height)
}

# * main
nmfDir <- args$nmfDir
files <- matchFiles(nmfDir = nmfDir, suffix = "sta.txt$")
nmfPrefix <- getNMFPrefix(files = files, extrastr = ".sta.txt$")
modules <- getNMFModules(nmfPrefix = nmfPrefix, files = NULL)
boxStaList <- readBoxSta(nmfDir = nmfDir, nmfPrefix = nmfPrefix,
  prefix2 = "n0")
statxt <- readStatxt(nmfDir = nmfDir,
  nmfPrefix = nmfPrefix, samples = "nmfPmat")
statH <- readStatH(nmfDir = nmfDir, nmfPrefix = nmfPrefix,
  prefix2 = "n0", samples = "nmfPmat")

# * plot box plot
withr::with_pdf(
  new = file.path(nmfDir, paste("nmf", args$tag, "statH.sta.pdf", sep = ".")),
  code = {
    p1 <- ggplot(statH, aes(ranks, contributes)) +
      geom_boxplot(aes(group = cut_width(ranks, 1), colour = ranks)) +
      theme_bw()
    print(p1)
    p2 <- ggplot(statH, aes(ranks, sparseness)) +
      geom_boxplot(aes(group = cut_width(ranks, 1), colour = ranks)) +
      theme_bw()
    print(p2)
    p3 <- ggplot(statH, aes(ranks, entropy)) +
      geom_boxplot(aes(group = cut_width(ranks, 1), colour = ranks)) +
      theme_bw()
    print(p3)
  }
)

# * plot entr 2 sparse
measure <- "Median"
entr_stats <- boxStaList$entropy
sparse_stats <- boxStaList$sparseness
module_nms <- paste("Module", colnames(entr_stats))

entr_values <- entr_stats[measure, ]
sparse_values <- sparse_stats[measure, ]
df_plot <- rbind(entr_values, sparse_values)
rownames(df_plot) <- c("Entropy", "Sparseness")
colnames(df_plot) <- module_nms

# * order the columns based on the size of module
## order.mod <- order(as.numeric(entr_stats), decreasing = FALSE)
## df_plot <- df_plot[, order.mod]

## # * plot 
## mat_plot <- as.matrix(t(df_plot))
## cols <- c("blue", "red")
## ## TODO: update to show two lines using two y-axis.
## withr::with_pdf(
##   new = file.path(nmfDir,
##     paste("nmf", args$tag, "entropy_sparse", "pdf",sep = ".")),
##   code = {
##     matplot(mat_plot, type = c("b"), pch = 19, col = cols,  xaxt = "n",
##       xlab = paste("Number of Modules for", args$tag),
##       ylab = "Measurements",
##       main = paste("NMF for", args$tag))
##     axis(side = 1, at = seq_len(nrow(mat_plot)), labels = rownames(mat_plot))
##     legend("topright", legend = c("Entropy", "Sparseness"), col = cols, pch = 19)
##   }, width = 10, height = 7)


