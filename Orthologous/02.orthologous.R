#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggradar"))

packdir <- "./package/R"
import::from(.from = "colors.R", .directory = packdir, ArchRPalettes)

# * function
posX_adjust <- function(x) {
    normalize_to_0_1 = (x - min(x)) / (max(x) - min(x))
    ratio_to_0_1 = x / sum(x)
    min_val <- pmin(normalize_to_0_1, ratio_to_0_1)
    return(min_val)
}
concatAllenAnnot1 <- function(allenAnnot, space2="_", slash2="-") {
  r <- gsub("/", slash2, allenAnnot)
  r <- gsub(" ", space2, r)
  r <- gsub("-", ".", r)
  return(r)
}

# * load data
peak.anno <- readRDS("./CEMBA_data/rds/update_pmat/peak.anno.rds")
peak.orthos <- readRDS("./CEMBA_data/rds/update_pmat/peak.orthos.rds")
peak.CAcons <- readRDS("./CEMBA_data/rds/update_pmat/peak.CAcons.rds")
peak.CAdiver <- readRDS("./CEMBA_data/rds/update_pmat/peak.CAdiver.rds")
peak.subclass <- readRDS("./CEMBA_data/rds/update_pmat/peak.subclass.rds")
celltype.dict <- readRDS("./CEMBA_data/rds/update_pmat/celltype.dict.rds")
L2.order <- readRDS("./CEMBA_data/rds/update_pmat/L2.order.rds")

# * mouse specific and orthologous
peak.anno %>% filter(!cCRE %in% peak.orthos$cCRE) -> peak.anno.specific
peak.anno %>% filter(cCRE %in% peak.orthos$cCRE) -> peak.anno.orthos
peak.anno %>% filter(cCRE %in% peak.CAcons$cCRE) -> peak.anno.CAcons
peak.anno %>% filter(cCRE %in% peak.CAdiver$cCRE) -> peak.anno.CAdiver

# * integrate data 
peak.anno.specific %>% group_by(AnnoSuperFam1) %>% dplyr::summarise(count=n()) %>% mutate(ratio=count/sum(count)) %>% 
        mutate(label=paste(AnnoSuperFam1,scales::percent(as.numeric(ratio), 0.1),sep=', ')) %>% 
        mutate(Ratio=str_replace(label,", 0.0%", "")) %>% column_to_rownames("AnnoSuperFam1") -> peak.specific.AnnoSuper
peak.anno.orthos %>% group_by(AnnoSuperFam1) %>% dplyr::summarise(count=n()) %>% mutate(ratio=count/sum(count)) %>% 
        mutate(label=paste(AnnoSuperFam1,scales::percent(as.numeric(ratio), 0.1),sep=', ')) %>% 
        mutate(Ratio=str_replace(label,", 0.0%", "")) %>% column_to_rownames("AnnoSuperFam1") -> peak.orthos.AnnoSuper
peak.anno.CAcons %>% group_by(AnnoSuperFam1) %>% dplyr::summarise(count=n()) %>% mutate(ratio=count/sum(count)) %>% 
        mutate(label=paste(AnnoSuperFam1,scales::percent(as.numeric(ratio), 0.1),sep=', ')) %>% 
        mutate(Ratio=str_replace(label,", 0.0%", "")) %>% column_to_rownames("AnnoSuperFam1") -> peak.CAcons.AnnoSuper
peak.anno.CAdiver %>% group_by(AnnoSuperFam1) %>% dplyr::summarise(count=n()) %>% mutate(ratio=count/sum(count)) %>% 
        mutate(label=paste(AnnoSuperFam1,scales::percent(as.numeric(ratio), 0.1),sep=', ')) %>% 
        mutate(Ratio=str_replace(label,", 0.0%", "")) %>% column_to_rownames("AnnoSuperFam1") -> peak.CAdiver.AnnoSuper

if(identical(rownames(peak.specific.AnnoSuper), rownames(peak.orthos.AnnoSuper))) {
    Anno_specVSortho_count <- rbind(peak.specific.AnnoSuper$count, peak.orthos.AnnoSuper$count)
    rownames(Anno_specVSortho_count) <- c("Mouse Spec.", "Orthologous")
    colnames(Anno_specVSortho_count) <- rownames(peak.specific.AnnoSuper)
    Anno_specVSortho_count %>% as.data.frame %>% mutate(cCREsNumber=rowSums(across(where(is.numeric)))) %>% rownames_to_column() -> Anno_sepcVSortho_totalcount
    Anno_specVSortho_ratio <- rbind(peak.specific.AnnoSuper$ratio, peak.orthos.AnnoSuper$ratio)
    rownames(Anno_specVSortho_ratio) <- c("Mouse Spec.", "Orthologous")
    colnames(Anno_specVSortho_ratio) <- rownames(peak.specific.AnnoSuper)
    Anno_specVSortho_ratio %>% as.data.frame %>% rownames_to_column("Set") -> Anno_specVSortho_ratio
    rownames(Anno_specVSortho_ratio) <- Anno_specVSortho_ratio$Set
}

if(identical(rownames(peak.CAcons.AnnoSuper), rownames(peak.CAdiver.AnnoSuper))){
    Anno_consVSdiver_count <- rbind(peak.CAcons.AnnoSuper$count, peak.CAdiver.AnnoSuper$count)
    rownames(Anno_consVSdiver_count) <- c("CAcons", "CAdiver")
    colnames(Anno_consVSdiver_count) <- rownames(peak.CAcons.AnnoSuper)
    Anno_consVSdiver_count %>% as.data.frame %>% mutate(cCREsNumber=rowSums(across(where(is.numeric)))) %>% rownames_to_column() -> Anno_consVSdiver_totalcount
    Anno_consVSdiver_ratio <- rbind(peak.CAcons.AnnoSuper$ratio, peak.CAdiver.AnnoSuper$ratio)
    rownames(Anno_consVSdiver_ratio) <- c("CAcons", "CAdiver")
    colnames(Anno_consVSdiver_ratio) <- rownames(peak.orthos.AnnoSuper)
    Anno_consVSdiver_ratio %>% as.data.frame %>% rownames_to_column("Set") -> Anno_consVSdiver_ratio
    rownames(Anno_consVSdiver_ratio) <- Anno_consVSdiver_ratio$Set
}

# * ortho. vs spec. barplot
cols_radar <- c("#46B8DA", "#EEA236")
names(cols_radar) <- c("Mouse Spec.", "Orthologous")
pt_bar_spec <- ggplot(Anno_sepcVSortho_totalcount, aes(x=as.factor(rowname), y=cCREsNumber)) + geom_bar(stat="identity", width=0.6, color=cols_radar[Anno_sepcVSortho_totalcount$rowname], 
    fill=cols_radar[Anno_sepcVSortho_totalcount$rowname]) + theme_classic() + theme(legend.position="none") + xlab("") + ylab("Number of cCREs")
ggasve(file="plot/Fig5a.specific.barplot.pdf", pt_bar_spec, width=8, height=8)

# * spec. pie plot
peak.specific.AnnoSuper$name <- rownames(peak.specific.AnnoSuper)
# Get the positions
peak.specific.AnnoSuper %>% mutate(csum=rev(cumsum(rev(count))), pos=count/2+lead(csum, 1), pos=if_else(is.na(pos), count/2, pos)) %>% 
    mutate(adjx=pmin(abs((pos-max(csum)*0)/max(csum)), abs((pos-max(csum)*0.25)/max(csum)), abs((pos-max(csum)*0.5)/max(csum)),
        abs((pos-max(csum)*0.75)/max(csum)), abs((pos-max(csum)*1)/max(csum)))/ratio) %>% mutate(adjx1=posX_adjust(adjx)) -> peak.specific.AnnoSuper
peak.specific.AnnoSuper %>% filter(pos > max(csum)/2) -> peak.specific.AnnoSuper.left
peak.specific.AnnoSuper %>% filter(pos <= max(csum)/2) -> peak.specific.AnnoSuper.right
# set colors
colors_specific.AnnoSuper <- ArchRPalettes[["stallion"]][1:nrow(peak.specific.AnnoSuper)]
names(colors_specific.AnnoSuper) <- rownames(peak.specific.AnnoSuper)
# plot pie
start_pos_AnnoSuper <- pi/9
adj.factor.right_AnnoSuper <- rep(1, nrow(peak.specific.AnnoSuper.right))
#adj.factor.right_AnnoSuper <- c(4,3,-1.5,3,-2,4,3,1.5)
adj.factor.left_AnnoSuper <- rep(1, nrow(peak.specific.AnnoSuper.left))
#adj.factor.left_AnnoSuper <- c(4,0.5,0.5,0.5,1)
pt_pie_specAnnoSuper <- ggplot(peak.specific.AnnoSuper, aes(x="" , y=count, fill=fct_inorder(name))) + geom_col(width=1, color=1, linewidth=0.1) + 
        coord_polar(theta="y", start=start_pos_AnnoSuper) + scale_fill_manual(values=colors_specific.AnnoSuper) + 
        geom_text_repel(data=peak.specific.AnnoSuper.left, aes(x=1.5, y=pos, label=Ratio), size=4.5, nudge_x=peak.specific.AnnoSuper.left$adjx1+start_pos_AnnoSuper/adj.factor.left_AnnoSuper, show.legend=FALSE, 
            segment.curvature=-1e-20, segment.angle=30, box.padding=0.5, min.segment.length=0.01, direction="x", hjust=0.5, segment.ncp=3) + 
        geom_text_repel(data=peak.specific.AnnoSuper.right, aes(x=1.5, y=pos, label=Ratio), size=4.5, nudge_x=peak.specific.AnnoSuper.right$adjx1+start_pos_AnnoSuper/adj.factor.right_AnnoSuper, show.legend=FALSE, 
                segment.curvature=-1e-20, segment.angle=30, box.padding=0.5, min.segment.length=0.01, direction="y", hjust=0, segment.ncp=3) + 
        guides(fill="none") + theme_void()
ggsave(file="./plot/FigS14a.specific.pie.pdf", pt_pie_specAnnoSuper, width=8, height=8)

# * ortho. vs spec. radar plot
Anno_specVSortho_count %>% as.data.frame %>% rownames_to_column("Set") -> Anno_specVSortho_count
rownames(Anno_specVSortho_count) <- Anno_specVSortho_count$Set
pt_radar_SO <- ggradar(Anno_specVSortho_ratio, grid.min=0.0, grid.mid=0.45/2, grid.max=0.45, 
    background.circle.colour="white",
    axis.line.colour="gray60",
    gridline.mid.colour="gray60",
    group.colours=cols_radar,
    values.radar=c("0%", "", "45%"),
    legend.title="Set",
    legend.position="bottom",
    grid.label.size=5, 
    axis.label.size=5,
    group.point.size=3)
ggsave(file="plot/Fig5b.radar.pdf", pt_radar_SO, width=8, height=8)

# * CA cons. vs diver. barplot
cols_radar <- c("#46B8DA", "#EEA236")
names(cols_radar) <- c("CAcons", "CAdiver")
pt_bar_ortho <- ggplot(Anno_consVSdiver_totalcount, aes(x=as.factor(rowname), y=cCREsNumber)) + geom_bar(stat="identity", width=0.6, color=cols_radar, 
    fill=cols_radar[c("CAcons", "CAdiver")]) + theme_classic() + theme(legend.position="none") + xlab("") + ylab("Number of cCREs")
ggsave(pt_bar_ortho, file="./plot/FigS13a.ortho.barplot.pdf", width=8, height=8)

# * radar plot
pt_radar_CD <- ggradar(Anno_consVSdiver_ratio, grid.min=0.0, grid.mid=0.45/2, grid.max=0.45, 
    background.circle.colour="white",
    axis.line.colour="gray60",
    #gridline.min.colour="gray60",
    gridline.mid.colour="gray60",
    #gridline.max.colour="gray60",
    #gridline.min.linetype = 1,
    #gridline.mid.linetype = 1,
    #gridline.max.linetype = 1,
    group.colours=cols_radar,
    values.radar=c("0%", "", "45%"),
    legend.title="Set",
    legend.position="bottom",
    grid.label.size=5, 
    axis.label.size=5,
    group.point.size=3)
ggsave(file="plot/FigS13c.ortho.radar.pdf", pt_radar_CD, width=8, height=8)

# * across subclass
peak.subclass.CAcons <- peak.subclass[,which(colnames(peak.subclass) %in% peak.anno.CAcons$cCRE)]
peak.subclass.CAdiver <- peak.subclass[,which(colnames(peak.subclass) %in% peak.anno.CAdiver$cCRE)]
peak.subclass.CAcons.sumPeak <- apply(peak.subclass.CAcons, 1, sum)
peak.subclass.CAcons.sumPeak <- peak.subclass.CAcons.sumPeak[match(L2.order$subclass, names(peak.subclass.CAcons.sumPeak))]
peak.subclass.CAcons.sumSubclass <- apply(peak.subclass.CAcons, 2, sum)
peak.subclass.CAdiver.sumPeak <- apply(peak.subclass.CAdiver, 1, sum)
peak.subclass.CAdiver.sumPeak <- peak.subclass.CAdiver.sumPeak[match(L2.order$subclass, names(peak.subclass.CAdiver.sumPeak))]
peak.subclass.CAdiver.sumSubclass <- apply(peak.subclass.CAdiver, 2, sum)

if(identical(names(peak.subclass.CAcons.sumPeak), names(peak.subclass.CAdiver.sumPeak))){
    sumPeak_CAcons_CAdiver <- data.frame(
        CAcons = peak.subclass.CAcons.sumPeak,
        CAdiver = peak.subclass.CAdiver.sumPeak,
        CAdiverNeg = -peak.subclass.CAdiver.sumPeak
    )   
}
sumPeak_CAcons_CAdiver$mainclass <- celltype.dict$mainclass[match(rownames(sumPeak_CAcons_CAdiver), celltype.dict$AllenAnnotConcat)]

pt_hist_CAconsVSdiver <- ggplot(sumPeak_CAcons_CAdiver, aes(x=x)) + 
    geom_histogram(aes(x=CAcons, y=after_stat(count)), fill="#69b3a2") + annotate("text", label="CA cons.", x=90000, y=20, size=6, colour="#69b3a2") + 
    geom_histogram(aes(x=CAdiver, y=-after_stat(count)), fill= "#404080") + annotate("text", label="CA diver.", x=90000, y=-20, size=6, colour="#404080")+ 
    ylab("Number of cell subclasses") + xlab("Number of CAcons or CAdiver cCREs in each subclass") + theme_classic()

ggsave(pt_hist_CAconsVSdiver, file="./plot/FigS13d.CAconsVSdiver.hist.pdf", width=8, height=8)

peak.subclass.CAcons.sumSubclass %>% as.data.frame %>% rownames_to_column %>% rename(Subclasses=".") -> peak.subclass.CAcons.sumSubclass
pt_hist_CAcons.cellspec <- ggplot(peak.subclass.CAcons.sumSubclass, aes(x=Subclasses)) + geom_histogram(fill="#69b3a2", color="#e9ecef") + 
    theme_classic() + xlab("Number of cell subclasses captured a same cCRE") + ylab("Number of cCREs") +
    annotate("text", label="CA cons.", x=100, y=40000, size=6, colour="black")

peak.subclass.CAdiver.sumSubclass %>% as.data.frame %>% rownames_to_column %>% rename(Subclasses=".") -> peak.subclass.CAdiver.sumSubclass
pt_hist_CAdiver.cellspec <- ggplot(peak.subclass.CAdiver.sumSubclass, aes(x=Subclasses)) + geom_histogram(fill="#69b3a2", color="#e9ecef") + 
    theme_classic() + xlab("Number of cell subclasses captured a same cCRE") + ylab("Number of cCREs") +
    annotate("text", label="CA diver.", x=100, y=40000, size=6, colour="black")
ggsave(pt_hist_CAcons.cellspec, file="./plot/FigS13e.CAcons.cellspec.hist.pdf", width=8, height=8)
ggsave(pt_hist_CAdiver.cellspec, file="./plot/FigS13f.CAdiver.cellspec.hist.pdf", width=8, height=8)

# * specific vs CAcons vs CAdiver across subclass
peak.subclass.data <- t(peak.subclass)
peak.subclass.data %>% as.data.frame %>% rownames_to_column("cCRE") -> peak.subclass.data
rownames(peak.subclass.data) <- peak.subclass.data$cCRE
peak.subclass.data %>% mutate(Peak_class=case_when(cCRE %in% peak.anno.specific$cCRE ~ "Mouse spec.", cCRE %in% peak.anno.CAcons$cCRE ~ "CA cons.", cCRE %in% peak.anno.CAdiver$cCRE ~ "CA diver.")) -> peak.subclass.data
peak.subclass.data %>% gather(Subclass, ifpeak, rownames(peak.subclass)) -> peak.subclass.data.long
peak.subclass.data.long %>% group_by(Peak_class, Subclass) %>% dplyr::summarise(count=sum(ifpeak)) -> peak.subclass.data.long.count
peak.subclass.data.long.count %>% group_by(Subclass) %>% mutate(sumpeak=sum(count)) -> peak.subclass.data.long.sumpeak
peak.subclass.data.long.sumpeak %>% mutate(fraction=count/sumpeak) -> peak.subclass.data.long.fraction
peak.subclass.data.long.fraction$mainclass <- celltype.dict$mainclass[match(peak.subclass.data.long.fraction$Subclass, celltype.dict$AllenAnnotConcat)]
peak.subclass.data.long_cons.diver <- peak.subclass.data.long %>% filter(Peak_class != "Mouse spec.")
peak.subclass.data.long_cons.diver %>% group_by(Peak_class, Subclass) %>% dplyr::summarise(count=sum(ifpeak)) -> peak.subclass.data.long.count_cons.diver
peak.subclass.data.long.count_cons.diver %>% group_by(Subclass) %>% mutate(sumpeak=sum(count)) -> peak.subclass.data.long.sumpeak_cons.diver
peak.subclass.data.long.sumpeak_cons.diver %>% mutate(fraction=count/sumpeak) -> peak.subclass.data.long.fraction_cons.diver
peak.subclass.data.long_cons.diver_totalcount <- peak.subclass.data.long.fraction %>% filter(Peak_class != "Mouse spec.")
if(identical(peak.subclass.data.long.fraction_cons.diver$count, peak.subclass.data.long_cons.diver_totalcount$count)){
    peak.subclass.data.long.fraction_cons.diver$totalpeak <- peak.subclass.data.long_cons.diver_totalcount$sumpeak
}
peak.subclass.data.long.fraction_cons <- peak.subclass.data.long.fraction_cons.diver %>% filter(Peak_class == "CA cons.") %>% arrange(desc(fraction))
peak.subclass.data.long.fraction_cons.diver$Subclass <- factor(peak.subclass.data.long.fraction_cons.diver$Subclass, levels=peak.subclass.data.long.fraction_cons$Subclass)

# cut subclasses to 5 groups based on total peaks
peak.subclass.data.long.fraction_cons.diver %>% filter(totalpeak<=quantile(peak.subclass.data.long.fraction_cons.diver$totalpeak,0.2)) -> peak.subclass.data.long.fraction_cons.diver_group1
peak.subclass.data.long.fraction_cons.diver %>% filter(totalpeak>quantile(peak.subclass.data.long.fraction_cons.diver$totalpeak,0.2)&totalpeak<=quantile(peak.subclass.data.long.fraction_cons.diver$totalpeak,0.4)) -> peak.subclass.data.long.fraction_cons.diver_group2
peak.subclass.data.long.fraction_cons.diver %>% filter(totalpeak>quantile(peak.subclass.data.long.fraction_cons.diver$totalpeak,0.4)&totalpeak<=quantile(peak.subclass.data.long.fraction_cons.diver$totalpeak,0.6)) -> peak.subclass.data.long.fraction_cons.diver_group3
peak.subclass.data.long.fraction_cons.diver %>% filter(totalpeak>quantile(peak.subclass.data.long.fraction_cons.diver$totalpeak,0.6)&totalpeak<=quantile(peak.subclass.data.long.fraction_cons.diver$totalpeak,0.8)) -> peak.subclass.data.long.fraction_cons.diver_group4
peak.subclass.data.long.fraction_cons.diver %>% filter(totalpeak>quantile(peak.subclass.data.long.fraction_cons.diver$totalpeak,0.8)) -> peak.subclass.data.long.fraction_cons.diver_group5

cols_radar <- c("#46B8DA", "#EEA236")
names(cols_radar) <- c("CA cons.", "CA diver.")

plot_bar_consVSdiver <- function(data, group="group 1") {
    pt_bar_consVSdiver_group <- ggplot(data, aes(x=Subclass, y=fraction, fill=Peak_class)) + geom_bar(position="stack", stat="identity")
    pt_bar_consVSdiver_group_stack <- pt_bar_consVSdiver_group + scale_fill_manual("Class", values=cols_radar) + 
        scale_y_continuous(expand = expansion(0)) + 
        ylab("Fraction of cCREs class in orthologous cCRE of each subclass") + xlab(paste("Subclass in total peak number quantile", group, sep="\n")) + theme_bw() + 
        theme(legend.position="bottom", axis.text.x=element_text(colour="black", size=10, angle=90, vjust=0.5, hjust=1), 
            axis.text.y=element_text(colour="black", size=10), axis.title=element_text(colour="black", size=16)) 
    return(pt_bar_consVSdiver_group_stack)    
}
pt_bar_consVSdiver_group1_stack <- plot_bar_consVSdiver(peak.subclass.data.long.fraction_cons.diver_group1, group="group 1")
pt_bar_consVSdiver_group2_stack <- plot_bar_consVSdiver(peak.subclass.data.long.fraction_cons.diver_group2, group="group 2")
pt_bar_consVSdiver_group3_stack <- plot_bar_consVSdiver(peak.subclass.data.long.fraction_cons.diver_group3, group="group 3")
pt_bar_consVSdiver_group4_stack <- plot_bar_consVSdiver(peak.subclass.data.long.fraction_cons.diver_group4, group="group 4")
pt_bar_consVSdiver_group5_stack <- plot_bar_consVSdiver(peak.subclass.data.long.fraction_cons.diver_group5, group="group 5")

pt_bar_consVSdiver_groupall <- ggarrange(pt_bar_consVSdiver_group1_stack, pt_bar_consVSdiver_group2_stack, pt_bar_consVSdiver_group3_stack, pt_bar_consVSdiver_group4_stack, pt_bar_consVSdiver_group5_stack, 
    labels=c("A", "B", "C", "D", "E"), ncol=5, nrow=1, align="h", common.legend = TRUE, legend = "bottom")
ggsave(pt_bar_consVSdiver_groupall, file="./plot/figS13b.consVSdiver_bar.pdf", width=24, height=10)




