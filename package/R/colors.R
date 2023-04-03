#' The color palettes from SnapATAC
#' @export
SnapATACPalette <- c("grey", "#E31A1C", "#FFD700", "#771122",
  "#777711", "#1F78B4", "#68228B", "#AAAA44",
  "#60CC52", "#771155", "#DDDD77", "#774411",
  "#AA7744", "#AA4455", "#117744", "#000080",
  "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9",
  "#BC80BD", "#FFED6F", "#7FC97F", "#BEAED4",
  "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
  "#BF5B17", "#666666", "#1B9E77", "#D95F02",
  "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A",
  "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
  "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928",
  "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4",
  "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
  "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8",
  "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
  "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A",
  "#984EA3", "#FFFF33", "#A65628", "#F781BF",
  "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA",
  "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
  "#FCCDE5")

#' List of color palettes that can be used in plots
#'
#' A collection of some original and some borrowed color palettes to
#' provide appealing color aesthetics for plots in ArchR.
#'
#' Comments in ArchR:
#' DISCLOSURE: This is a collection of palettes that includes some
#' original palettes and some palettes originally
#' implemented by others in other packages.
#' They are included here for convenience because they help improve
#' plot aesthetics.
#' NOTE: all palettes included in the "Primarily Continuous Palettes"
#' section should also work for discrete usage but not vice versa.
#' Each continuous palette has been ordered by color to generate a
#' visually appealing discrete palette.
#'
#' @export
ArchRPalettes <- list(

  #---------------------------------------------------------------
  # Primarily Discrete Palettes
  #---------------------------------------------------------------

  # 20-colors
  stallion = c("1" = "#D51F26", "2" = "#272E6A", "3" = "#208A42", "4"
    = "#89288F", "5" = "#F47D2B", "6" = "#FEE500", "7" = "#8A9FD1", "8"
    = "#C06CAB", "19" = "#E6C2DC",
    "10" = "#90D5E4", "11" = "#89C75F", "12" = "#F37B7D", "13" =
      "#9983BD", "14" = "#D24B27", "15" = "#3BBCA8", "16" = "#6E4B9E",
    "17" = "#0C727C", "18" = "#7E1416", "9" = "#D8A767", "20" =
      "#3D3D3D"),

  stallion2 = c("1" = "#D51F26", "2" = "#272E6A", "3" = "#208A42", "4"
    = "#89288F", "5" = "#F47D2B", "6" = "#FEE500", "7" = "#8A9FD1",
    "8" = "#C06CAB", "19" = "#E6C2DC",
    "10" = "#90D5E4", "11" = "#89C75F", "12" = "#F37B7D", "13" =
      "#9983BD", "14" = "#D24B27", "15" = "#3BBCA8", "16" = "#6E4B9E",
    "17" = "#0C727C", "18" = "#7E1416", "9" = "#D8A767"),

  calm = c("1" = "#7DD06F", "2" = "#844081", "3" = "#688EC1", "4" =
    "#C17E73", "5" = "#484125", "6" = "#6CD3A7", "7" = "#597873", "8"
  = "#7B6FD0", "9" = "#CF4A31", "10" = "#D0CD47",
  "11" = "#722A2D", "12" = "#CBC594", "13" = "#D19EC4", "14" =
    "#5A7E36", "15" = "#D4477D", "16" = "#403552", "17" = "#76D73C",
  "18" = "#96CED5", "19" = "#CE54D1", "20" = "#C48736"),

  kelly = c("1" = "#FFB300", "2" = "#803E75", "3" = "#FF6800", "4" =
    "#A6BDD7", "5" = "#C10020", "6" = "#CEA262", "7" = "#817066", "8"
  = "#007D34", "9" = "#F6768E", "10" = "#00538A",
  "11" = "#FF7A5C", "12" = "#53377A", "13" = "#FF8E00", "14" =
    "#B32851", "15" = "#F4C800", "16" = "#7F180D", "17" = "#93AA00",
  "18" = "#593315", "19" = "#F13A13", "20" = "#232C16"),

  # 16-colors
  bear = c("1" = "#faa818", "2" = "#41a30d", "3" = "#fbdf72", "4" =
    "#367d7d", "5" = "#d33502", "6" = "#6ebcbc", "7" = "#37526d",
  "8" = "#916848", "9" = "#f5b390", "10" = "#342739", "11" =
    "#bed678", "12" = "#a6d9ee", "13" = "#0d74b6",
  "14" = "#60824f", "15" = "#725ca5", "16" = "#e0598b"),

  # 15-colors
  ironMan = c("9" = '#371377', "3" = '#7700FF', "2" = '#9E0142', "10"
    = '#FF0080', "14" = '#DC494C', "12" = "#F88D51", "1" = "#FAD510",
    "8" = "#FFFF5F", "4" = '#88CFA4',
    "13" = '#238B45', "5" = "#02401B", "7" = "#0AD7D3", "11" =
      "#046C9A", "6" = "#A2A475", "15" = 'grey35'),

  circus = c("1" = "#D52126", "2" = "#88CCEE", "3" = "#FEE52C", "4" =
    "#117733", "5" = "#CC61B0", "6" = "#99C945", "7" = "#2F8AC4", "8"
  = "#332288",
  "9" = "#E68316", "10" = "#661101", "11" = "#F97B72", "12" =
    "#DDCC77", "13" = "#11A579", "14" = "#89288F", "15" = "#E73F74"),

  # 12-colors
  paired = c("9" = "#A6CDE2", "1" = "#1E78B4", "3" = "#74C476", "12" =
    "#34A047", "11" = "#F59899", "2" = "#E11E26",
  "10" = "#FCBF6E", "4" = "#F47E1F", "5" = "#CAB2D6", "8" =
    "#6A3E98", "6" = "#FAF39B", "7" = "#B15928"),

  # 11-colors
  grove = c("11" = "#1a1334", "9" = "#01545a", "1" = "#017351", "6" =
    "#03c383", "8" = "#aad962", "2" = "#fbbf45", "10" = "#ef6a32", "3" =
    "#ed0345", "7" = "#a12a5e", "5" = "#710162", "4" = "#3B9AB2"),

  # 7-colors
  summerNight = c("1" = "#2a7185", "2" = "#a64027", "3" = "#fbdf72",
    "4" = "#60824f", "5" = "#9cdff0", "6" = "#022336", "7" = "#725ca5"),

  # 5-colors
  zissou = c("1" = "#3B9AB2", "4" = "#78B7C5", "3" = "#EBCC2A", "5" =
    "#E1AF00", "2" = "#F21A00"), # wesanderson
  darjeeling = c("1" = "#FF0000", "2" = "#00A08A", "3" = "#F2AD00",
    "4" = "#F98400", "5" = "#5BBCD6"), # wesanderson
  rushmore = c("1" = "#E1BD6D", "5" = "#EABE94", "2" = "#0B775E", "4"
    = "#35274A", "3" = "#F2300F"), # wesanderson
  captain = c("1" = "grey", "2" = "#A1CDE1", "3" = "#12477C", "4" =
    "#EC9274", "5" = "#67001E"),

  #---------------------------------------------------------------
  # Primarily Continuous Palettes
  #---------------------------------------------------------------

  # 10-colors
  horizon = c("1" = "#000075", "4" = "#2E00FF", "6" = "#9408F7", "10"
    = "#C729D6", "8" = "#FA4AB5", "3" = "#FF6A95", "7" = "#FF8B74", "5"
    = "#FFAC53", "9" = "#FFCD32", "2" = "#FFFF60"),

  # 9-colors
  horizonExtra = c("1" = "#000436", "4" = "#021EA9", "6" = "#1632FB",
    "8" = "#6E34FC", "3" = "#C732D5", "9" = "#FD619D", "7" = "#FF9965",
    "5" = "#FFD32B", "2" = "#FFFC5A"),
  blueYellow = c("1" = "#352A86", "2" = "#343DAE", "3" = "#0262E0",
    "4" = "#1389D2", "5" = "#2DB7A3", "6" = "#A5BE6A", "7" = "#F8BA43",
    "8" = "#F6DA23", "9" = "#F8FA0D"),
  sambaNight = c("6" = "#1873CC", "2" = "#1798E5", "8" = "#00BFFF",
    "5" = "#4AC596", "1" = "#00CC00", "4" = "#A2E700", "9" = "#FFFF00",
    "7" = "#FFD200", "3" = "#FFA500"), # buencolors
  solarExtra = c("5" = "#3361A5", "7" = "#248AF3", "1" = "#14B3FF",
    "8" = "#88CEEF", "9" = "#C1D5DC", "4" = "#EAD397", "3" = "#FDB31A",
    "2" = "#E42A2A", "6" = "#A31D1D"), # buencolors
  whitePurple = c("9" = "#f7fcfd", "6" = "#e0ecf4", "8" = "#bfd3e6",
    "5" = "#9ebcda", "2" = "#8c96c6", "4" = "#8c6bb1", "7" = "#88419d",
    "3" = "#810f7c", "1" = "#4d004b"),
  whiteBlue = c("9" = "#fff7fb", "6" = "#ece7f2", "8" = "#d0d1e6", "5"
    = "#a6bddb", "2" = "#74a9cf", "4" = "#3690c0", "7" = "#0570b0", "3"
    = "#045a8d", "1" = "#023858"),
  whiteRed = c("1" = "white", "2" = "red"),
  comet = c("1" = "#E6E7E8", "2" = "#3A97FF", "3" = "#8816A7", "4" =
    "black"),

  # 7-colors
  greenBlue = c("4" = "#e0f3db", "7" = "#ccebc5", "2" = "#a8ddb5", "5"
    = "#4eb3d3", "3" = "#2b8cbe", "6" = "#0868ac", "1" = "#084081"),

  # 6-colors
  beach = c("4" = "#87D2DB", "1" = "#5BB1CB", "6" = "#4F66AF", "3" =
    "#F15F30", "5" = "#F7962E", "2" = "#FCEE2B"),

  # 5-colors
  coolwarm = c("1" = "#4858A7", "4" = "#788FC8", "5" = "#D6DAE1", "3"
    = "#F49B7C", "2" = "#B51F29"),
  fireworks = c("5" = "white", "2" = "#2488F0", "4" = "#7F3F98", "3" =
    "#E22929", "1" = "#FCB31A"),
  greyMagma = c("2" = "grey", "4" = "#FB8861FF", "5" = "#B63679FF",
    "3" = "#51127CFF", "1" = "#000004FF"),
  fireworks2 = c("5" = "black", "2" = "#2488F0", "4" = "#7F3F98", "3"
    = "#E22929", "1" = "#FCB31A"),
  purpleOrange = c("5" = "#581845", "2" = "#900C3F", "4" = "#C70039",
    "3" = "#FF5744", "1" = "#FFC30F")
)


#' @export
getClassColors <- function() {
  classColors <- c("#90D5E4", "#89C75F", "#F37B7D",
    "#D51F26", "#89288F", "#8A9FD1", "#C06CAB", "#D8A767")
  names(classColors) <- c("NN", "Glut", "Gaba",
    "Dopa", "Chol", "Gly", "Sero", "Nora")
  return(classColors)
}

#' @export
mainclassColorsV1 <- c(
  "GABA" = "#5B7C3D",
  "GLUT" = "#8C4E3C",
  "NN" = "#704E89"
)
#' @export
mainclassColorsV2 <- c(
  "GABA" = "#D51F26",
  "GLUT" = "#8DC360",
  "NN" = "#5BBCD6")

#' @export
classColors <- c(
  S4Vectors::Pairs("NN", "#90D5E4"),
  S4Vectors::Pairs("Glut", "#89C75F"),
  S4Vectors::Pairs("GABA", "#F37B7D"),
  S4Vectors::Pairs("Dopa", "#D51F26"),
  S4Vectors::Pairs("Nora", "#D8A767"),
  S4Vectors::Pairs("Sero", "#C06CAB"),
  S4Vectors::Pairs("GABA-Glut", "#8C4E3C"),
  S4Vectors::Pairs("GABA-Chol", "#89288F"),
  S4Vectors::Pairs("GABA-Dopa", "#5B7C3D"),
  S4Vectors::Pairs("GABA-Gly", "#8A9FD1"),
  S4Vectors::Pairs("GABA-Hist", "#000080"),
  S4Vectors::Pairs("Glut-Chol", "#704E89"),
  S4Vectors::Pairs("Glut-Dopa", "#FB8072")
)

#' @export
largeRegionColors <- c(
  S4Vectors::Pairs("Telencephalon", "#00688B"),
  S4Vectors::Pairs("Diencephalon", "#F15F30"),
  S4Vectors::Pairs("Midbrain", "#74E44B"),
  S4Vectors::Pairs("Hindbrain", "#788FC8"),
  S4Vectors::Pairs("Cerebellum", "#DEB34C"),
  S4Vectors::Pairs("Non-Telencephalon", "#999999")
)

#' Set under the help of getL2Colors, then
#' modify the light colors.
#' @export
L2Annot2Colors <- c(
  S4Vectors::Pairs("MGL", "#CD96CD"),
  S4Vectors::Pairs("VLMC", "#8B4789"),
  S4Vectors::Pairs("VPIA", "#8B1C62"),
  S4Vectors::Pairs("PER", "#8B668B"),
  S4Vectors::Pairs("VEC", "#EE30A7"),
  S4Vectors::Pairs("OPC", "#E22929"),
  S4Vectors::Pairs("IOL", "#C732D5"),
  S4Vectors::Pairs("OGC1", "#EE7AE9"),
  S4Vectors::Pairs("OGC2", "#CD69C9"),
  S4Vectors::Pairs("CHOR", "#8470FF"),
  S4Vectors::Pairs("EPTANN", "#8968CD"),
  S4Vectors::Pairs("BERG", "#7B68EE"),
  S4Vectors::Pairs("ASCTENT", "#FF34B3"),
  S4Vectors::Pairs("ASCNT", "#AB82FF"),
  S4Vectors::Pairs("RGL", "#9F79EE"),
  S4Vectors::Pairs("ASCTE1", "#FF83FA"),
  S4Vectors::Pairs("ASCTE2", "#CD2990"),
  S4Vectors::Pairs("ITRSPL23GL2", "#FB8861FF"),
  S4Vectors::Pairs("ITRSPL23GL1", "#EEA9B8"),
  S4Vectors::Pairs("RSPL4GL", "#CD919E"),
  S4Vectors::Pairs("ITL5GL", "#A0522D"),
  S4Vectors::Pairs("ITL45GL", "#8B4726"),
  S4Vectors::Pairs("ITL46GL", "#FF8C69"),
  S4Vectors::Pairs("NPGL", "#EE8262"),
  S4Vectors::Pairs("CTL6GL", "#EE7942"),
  S4Vectors::Pairs("L6bGL", "#FF6347"),
  S4Vectors::Pairs("CLAGL", "#8B4C39"),
  S4Vectors::Pairs("ETL5GL", "#CD7054"),
  ## S4Vectors::Pairs("PVGA2", "#009ACD"),
  S4Vectors::Pairs("PVGL", "#FC8D62"),
  S4Vectors::Pairs("ITL23GL", "#FF8247"),
  S4Vectors::Pairs("ITL56GL", "#CD6839"),
  S4Vectors::Pairs("HPFGL", "#8B636C"),
  S4Vectors::Pairs("CA1GL", "#8B5A2B"),
  S4Vectors::Pairs("CA3GL", "#D2B48C"),
  S4Vectors::Pairs("AMYGL", "#FA8072"),
  S4Vectors::Pairs("ITENTPIRGL2", "#8B3626"),
  S4Vectors::Pairs("ITENTPIRGL1", "#CD4F39"),
  S4Vectors::Pairs("RHPGL", "#CD853F"),
  S4Vectors::Pairs("HBGAGL", "#36648B"),
  S4Vectors::Pairs("HBGL", "#20B2AA"),
  S4Vectors::Pairs("CBGRGL", "#EE5C42"),
  S4Vectors::Pairs("CBGA", "#7A8B8B"),
  S4Vectors::Pairs("THGL", "#FFA54F"),
  S4Vectors::Pairs("PURKGA", "#046C9A"),
  S4Vectors::Pairs("MYGL", "#63B8FF"),
  S4Vectors::Pairs("MHGL", "#02401B"),
  S4Vectors::Pairs("STRGL", "#4F94CD"),
  S4Vectors::Pairs("MBDOP", "#352A86"),
  S4Vectors::Pairs("THGA", "#80B1D3"),
  S4Vectors::Pairs("ICGL", "#6C7B8B"),
  S4Vectors::Pairs("MBGA1", "#377EB8"),
  S4Vectors::Pairs("MSGA", "#458B74"),
  S4Vectors::Pairs("MBHBGLGA2", "#00538A"),
  S4Vectors::Pairs("MBGL", "#9AC0CD"),
  S4Vectors::Pairs("MBHBGLGA1", "#53868B"),
  S4Vectors::Pairs("MBGA2", "#9FB6CD"),
  S4Vectors::Pairs("MXDGA", "#1874CD"),
  S4Vectors::Pairs("D12MSN", "#104E8B"),
  S4Vectors::Pairs("STRPALGA", "#1C86EE"),
  S4Vectors::Pairs("PIRGA", "#0d74b6"),
  S4Vectors::Pairs("VIPGA", "#00688B"),
  S4Vectors::Pairs("LAMGA", "#1E90FF"),
  S4Vectors::Pairs("STRGA", "#66CDAA"),
  S4Vectors::Pairs("PVGA", "#00BFFF"),
  S4Vectors::Pairs("SSTGA", "#00B2EE"),
  S4Vectors::Pairs("LSXGA", "#0B775E"),
  S4Vectors::Pairs("AMYGAGL", "#2a7185"),
  S4Vectors::Pairs("HYGLGA", "#3B9AB2"),
  S4Vectors::Pairs("OBGA", "#01545a"),
  S4Vectors::Pairs("OBNBL", "#7AC5CD"),
  S4Vectors::Pairs("CRGL", "#68838B"),
  S4Vectors::Pairs("ITENTPIRGL3", "#EE9A49"),
  S4Vectors::Pairs("DGNBL", "#5CACEE"),
  S4Vectors::Pairs("LQ", "#BEBEBE")
)

#' @export
#' Ref:https://stackoverflow.com/questions/31574480/rgb-to-hex-converter
#' @rgbmat: 3 by n
rgb2hex <- function(rgbmat) {
  # function to apply to each column of input rgbmat
  ProcessColumn <- function(col) {
    rgb(rgbmat[1, col],
      rgbmat[2, col],
      rgbmat[3, col],
      maxColorValue = 255)
  }
  # Apply the function
  sapply(seq_len(ncol(rgbmat)), ProcessColumn)
}

#' @export
getL2Colors <- function(cembav2env) {
  hcL2 <- readRDS(cembav2env$L2HCFile)
  L2.order <- hcL2$labels[hcL2$order]
  atacMeta <- readRDS(cembav2env$allmetaFile)
  L2.L2Annot <- unique(atacMeta[, c("L2", "L2Annot2")])
  rownames(L2.L2Annot) <- L2.L2Annot$L2

  mcL2 <- rep("NN", nrow(L2.L2Annot))
  mcL2[grepl("GABA", L2.L2Annot$L2)] <- "GABA"
  mcL2[grepl("GLUT", L2.L2Annot$L2)] <- "GLUT"
  names(mcL2) <- L2.L2Annot$L2Annot2
  mcL2 <- mcL2[L2.order]
  colorL2.NN <- c(
    paste0("orchid", 1:4), paste0("plum", 1:4),
    paste0("maroon", 1:4), paste0("mediumpurple", 1:4),
    "mediumslateblue", "lightslateblue"
  )
  names(colorL2.NN) <- names(mcL2)[mcL2 %in% "NN"]
  colorL2.GLUT <- c(
    paste0("sienna", 1:4), paste0("salmon", 1:4),
    paste0("tomato", 1:4), paste0("tan", 1:4),
    "tan", "sienna", "salmon", paste0("pink", 1:4)
  )
  names(colorL2.GLUT) <- names(mcL2)[mcL2 %in% "GLUT"]

  colorL2.GABA <- c(
    paste0("deepskyblue", 1:4), paste0("dodgerblue", 1:4),
    paste0("aquamarine", 1:4), paste0("cadetblue", 1:4),
    paste0("lightblue", 1:4), paste0("steelblue", 1:4),
    paste0("slategray", 1:4), paste0("lightcyan", 1:4),
    "lightseagreen"
  )
  names(colorL2.GABA) <- names(mcL2)[mcL2 %in% "GABA"]
  colorL2 <- c(colorL2.NN, colorL2.GLUT, colorL2.GABA)
  colorL2 <- colorL2[L2.order]
  colorL2["LQ"] <- "gray"
  # to rgb
  colorL2 <- grDevices::col2rgb(colorL2)
  # to hex
  colorL2 <- rgb2hex(colorL2)
  names(colorL2) <- c(L2.order, "LQ")
  r <- data.frame(L2Annot2 = names(colorL2), color = colorL2)
  return(colorL2)
}

#' @export
getSubclassColors.1 <- function(subclass,
                                mainclass,
                                dark.nn = "deepskyblue4",
                                light.nn = "darkslategray2",
                                dark.glut = "darkolivegreen",
                                light.glut = "darkseagreen1",
                                dark.gaba = "indianred4",
                                light.gaba = "lightpink1",
                                bias = 1,
                                alpha = FALSE,
                                space = "rgb",
                                interpolate = "linear") {
  mainclass <- toupper(mainclass)
  subclass.list <- list(
    nn = subclass[grep("NN|NONN", mainclass)],
    glut = subclass[grep("GLUT", mainclass)],
    gaba = subclass[grep("GABA", mainclass)]
  )
  color.border.list <- list(
    nn = c(light.nn, dark.nn),
    glut = c(light.glut, dark.glut),
    gaba = c(light.gaba, dark.gaba)
  )
  color.subclass.list <- lapply(c("nn", "glut", "gaba"),
    function(i) {
      r <- grDevices::colorRampPalette(
        colors = color.border.list[[i]],
        bias = bias,
        space = space,
        interpolate = interpolate,
        alpha = alpha)(
          length(subclass.list[[i]]))
      return(r)
    })
  color.subclass <- unlist(color.subclass.list)
  names(color.subclass) <- unlist(subclass.list)
  return(color.subclass)
}

#' Use R colors firstly: colors()
#' https://r-graph-gallery.com/ggplot2-color.html
#' @export
getAllRcolors <- function(file = "../meta/R.color.pdf") {
 withr::with_pdf(
  new = file,
  code = {
    scales::show_col(colors(), labels = TRUE, cex_label = 0.5)
  },
  width = 20,
  height = 20)
}
