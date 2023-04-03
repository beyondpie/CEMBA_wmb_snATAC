library(SnapATAC)
library(rjson)
library(optparse)
library(pureRUtils)

op <- list(
  make_option(c("--config"), type = "character", default = "config.json"),
  make_option(c("--group"), type = "character", default = "NA")
)

args <- parse_args(OptionParser(option_list = op))

pureRUtils::checkArgsExistOrStop(args)
pureRUtils::checkFileExistOrStop(args$config)

config <- rjson::fromJSON(file = args$config)

groupType <- config$groupType

ndpL2 <- config$ndpL2
ndpL3 <- config$ndpL3

matDir <- file.path(config$matDir, groupType)
pureRUtils::prepareOutdir(matDir)
matType <- config$mat
gencode <- config$gencode
## group <- args$group
if(config$splitGroup > 0) {
  message("Split group based on pattern ,")
  group <- as.vector(stringr::str_split(args$group, ",", simplify = T))
} else {
  message("Use raw group: ", args$group)
  group <- c(args$group)
}
requireNextLevel <- config$requireNextLevel

# * functions
getL2Groups <- function(annot, groupType = "L2") {
  if (groupType == "L2Id") {
    group <- paste(annot$L1Id, annot$L2Id, sep = ".")
  }
  if (groupType == "MajorType") {
    group <- paste(annot$L1Label, annot$L2Label, sep = ".")
  }
  if (groupType == "L2CoId") {
    group <- annot$L2CoId
  }
  ## this is now current usage
  if (groupType == "L2") {
    group <- annot$L2
  }
  if(groupType == "L2Extra") {
    group <- annot$L2Extra
  }
  return(group)
}

getL3GroupsForL2GroupType <- function(annot, groupType = "L2") {
  if (groupType == "L2Id") {
    group <- NULL
  }
  if (groupType == "MajorType") {
    group <- annot$L3Label
  }
  if (groupType == "L2CoId") {
    group <- paste(annot$L2CoId, annot$L3Id, sep = ".")
  }
  ## this is now current usage.
  if (groupType == "L2") {
    group <- annot$L3
  }
  if (groupType == "L2Extra") {
    group <- annot$L3Extra
  }
  if (length(group) < 1) {
    message("No L3 Group match L2 group ID of ", l2GroupID)
    message("Return NULL")
    return(NULL)
  }
  return(group)
}

# * load annotFile and get group information
## annot <- readRDS(config$annotFile)

# * load snap file
if (gencode == "vM16") {
  message("Read gencode.vM16 snap file.")
  snap <- readRDS(config$gencodevM16snapFile)
} else {
  message("Read gencode.vM23 snap file.")
  snap <- readRDS(config$gencodevM23snapFile)
}
snap <- pureRUtils::orderSnap(snap)
annot <- snap@metaData
if(groupType == "L1") {
  message("Group type is L1.")
  l1Groups <- annot$L1
  outDir <- file.path(matDir, group)
  pureRUtils::prepareOutdir(outDir)
  l1 <- group
  message("Prepare ", matType, " for group: ", l1)
  invisible(pureRUtils::addMatToSnap(snap = snap,
                                     cellGroup = l1Groups,
                                     group = l1,
                                     outDir = outDir,
                                     matType = matType,
                                     gencode = gencode,
                                     ndp = ndpL2,
                                     ncores = config$ncores,
                                     coverFile = TRUE,
                                     prefix = l1))
  
} else {
  message("Group type is: ", groupType)
  l2Groups <- getL2Groups(annot = annot, groupType = groupType)

  # * main
  if (group == "NA") {
    uL2Groups <- unique(l2Groups)
  } else {
    uL2Groups <- c(group)
  }

  invisible(lapply(uL2Groups, function(l2) {
    if (groupType == "L2Id") {
      outDir <- matDir
    } else {
      outDir <- file.path(matDir, l2)
    }
    pureRUtils::prepareOutdir(outDir)
    message("Prepare ", matType, " for group: ", l2)
    invisible(pureRUtils::addMatToSnap(snap = snap,
                                       cellGroup = l2Groups,
                                       group = l2,
                                       outDir = outDir,
                                       matType = matType,
                                       gencode = gencode,
                                       ndp = ndpL2,
                                       ncores = config$ncores,
                                       coverFile = TRUE,
                                       prefix = l2))
    if (requireNextLevel < 1) {
      message("No requirment for next level.")
      return(NULL)
    }
    ## generate L3-level mat for L2CoId and Majortype
    l3Groups <- getL3GroupsForL2GroupType(
      annot = snap@metaData,
      groupType = groupType)
    if (is.null(l3Groups)) {
      message("No l3 groups found, skip.")
      return(NULL)
    }
    ## l3 labels only for the l2 group.
    l2Groups <- getL2Groups(annot = snap@metaData, groupType = groupType)
    localL3Groups <- l3Groups[l2Groups %in% l2]
    uL3Groups <- unique(localL3Groups)
    if (length(uL3Groups) < 2) {
      message("Only one group in L3 level, skip it for", l2)
      return(NULL)
    }
    invisible(lapply(uL3Groups, function(g) {
      message(paste("Prepare L3 group", g,
                    "of L2", l2, "for", matType, sep = " "))
      pureRUtils::addMatToSnap(snap = snap,
                               cellGroup = l3Groups,
                               group = g,
                               outDir = outDir,
                               matType = matType,
                               gencode = gencode,
                               ndp = ndpL3,
                               ncores = config$ncores,
                               coverFile = TRUE,
                               prefix = g)
    }))
    return(NULL)
  }))
}
