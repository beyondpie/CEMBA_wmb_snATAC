#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggradar"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("mixtools"))

packdir <- "./package/R"
import::from(.from = "colors.R", .directory = packdir, ArchRPalettes)

# * function
posX_adjust <- function(x) {
    normalize_to_0_1 = (x - min(x)) / (max(x) - min(x))
    ratio_to_0_1 = x / sum(x)
    min_val <- pmin(normalize_to_0_1, ratio_to_0_1)
    return(min_val)
}


smScatter_plot <- function(TEs.meta, cv, mean, title="Title", xlab="x-axis", ylab="y-axis", alpha=0){
    buylrd <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
        "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
    g1 <- ggplot(data=TEs.meta, aes_string(cv, mean)) + stat_density2d(aes(fill=..density..^0.25), geom="tile", contour=FALSE, n=256) + 
        scale_fill_gradientn(colours=buylrd) 
    g2 <- g1 + geom_point(alpha=alpha, shape=20)
    g3 <- g2 + labs(title=title, x=xlab, y=ylab)
    g4 <- g3 + theme(plot.title=element_text(face="bold.italic", size="14", color="black", hjust=0.5), 
        axis.title=element_text(face="bold", size=10, color="black"),
        axis.text=element_text(face="plain", size=8, color="black"),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid.major=element_blank(),legend.position="none",
        panel.spacing=unit(1, "cm"),
        plot.margin=margin(1,1,1,1,"cm"))
    return(g4)
}

# * load data
Anno_sepcVSortho_totalcount <- readRDS("Anno_sepcVSortho_totalcount.rds")
peak.specific.AnnoSuper <- readRDS("peak.specific.AnnoSuper.rds")
Anno_specVSortho_count <- readRDS("Anno_specVSortho_count.rds")
Anno_consVSdiver_totalcount <- readRDS("Anno_consVSdiver_totalcount.rds")
Anno_consVSdiver_ratio <- readRDS("Anno_consVSdiver_ratio.rds")
Anno_specVSortho_ratio <- readRDS("Anno_specVSortho_ratio.rds")
sumPeak_CAcons_CAdiver <- readRDS("sumPeak_CAcons_CAdiver.rds")
peak.subclass.CAcons.sumSubclass <- readRDS("peak.subclass.CAcons.sumSubclass.rds")
peak.subclass.CAdiver.sumSubclass <- readRDS("peak.subclass.CAdiver.sumSubclass.rds")
peak.subclass.data.long.fraction_cons.diver <- readRDS("peak.subclass.data.long.fraction_cons.diver.rds")
normalmix <- readRDS("normalmix.rds")
celltype.dict <- readRDS("celltype.dict.rds")
L2.order <- readRDS("L2.order.rds")
peak.subclass.AnnoRatio <- readRDS("peak.subclass.AnnoRatio.rds")
peak.subclass.TEfrac <- readRDS("peak.subclass.TEfrac.rds")
TEs.count.TE.meta.sum <- readRDS("TEs.count.TE.meta.sum.rds")
TEs.count.invariableTE.meta <- readRDS("TEs.count.invariableTE.meta.rds")

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


# * subclass cCRE annotation heatmap
celltype.dict_order <- celltype.dict[match(colnames(peak.subclass.AnnoRatio), celltype.dict$AllenAnnotConcat), ]
ha_col <- HeatmapAnnotation(mainclass=celltype.dict_order[["mainclass"]], class=celltype.dict_order[["classv2"]], which = "column",
    col=list(mainclass=setNames(brewer.pal(name="Dark2", n=3), unique(celltype.dict_order[,"mainclass"])),
        class=setNames(brewer.pal(name="Paired", n=12), unique(celltype.dict_order[,"classv2"]))))
cols_order <- match(L2.order$subclass, colnames(peak.subclass.AnnoRatio))
col_fun = circlize::colorRamp2(c(0, 0.05, 0.5), c("blue", "white", "red"))
pt_AnnoSubclass_heat <- Heatmap(peak.subclass.AnnoRatio, name="Percentage", col=col_fun, 
    column_order=cols_order, row_order=c(1:13),
    show_row_names=TRUE, 
    row_names_gp=gpar(fontsize=8),
    show_column_names=TRUE,
    column_names_gp=gpar(fontsize=4),
    top_annotation=ha_col)
ggsave(pt_AnnoSubclass_heat, file="./plot/figS14d.Subclass.cCRE.heatmap.pdf", width=16, height=9)

# * TE fraction in each subclass histogram plot
cols_mainclass <- celltype.dict %>% select(mainclass, mainclassColorv2) %>% distinct %>% select(mainclassColorv2) %>% unlist
names(cols_mainclass) <- celltype.dict %>% select(mainclass, mainclassColorv2) %>% distinct %>% select(mainclass) %>% unlist

pt_hist_TEfrac <- ggplot(peak.subclass.TEfrac, aes(x=ratio, fill=mainclass, color=mainclass)) + geom_histogram(alpha=0.6, position="identity") +
    theme_classic() + scale_color_manual(values=cols_mainclass) + scale_fill_manual(values=cols_mainclass) + labs(x="Fraction of cCREs overlpped with TEs in subclass", y="# of cell subclasses") + 
    annotate("text", x=0.05, y=11, label="Cell subclass")

pt_hist_TEfrac_facets <- ggplot(peak.subclass.TEfrac, aes(x=ratio)) + geom_histogram(aes(color=mainclass, fill=mainclass)) + facet_grid(mainclass ~ .) +
    theme_classic() + scale_color_manual(values=cols_mainclass) + scale_fill_manual(values=cols_mainclass) + theme(legend.position="none") + 
    labs(x="Fraction of cCREs overlpped with TEs in subclass", y="# of cell subclasses")

# * highTE cell subclasses
pt_hist_TEfrac <- ggplot(peak.subclass.TEfrac, aes(x=ratio, fill=mainclass, color=mainclass)) + geom_histogram(alpha=0.6, position="stack", binwidth=(max(peak.subclass.TEfrac$ratio)-min(peak.subclass.TEfrac$ratio))/150) +
  theme_classic() +scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2") + labs(x="Fraction of open cCREs overlap with TEs in subclass", y="Number of cell subclasses") + 
  annotate("text", x=0.05, y=11, label="Cell subclass")
pt_hist_fit <- pt_hist_TEfrac + mapply(function(mean, sd, lambda, n, binwidth) {
  stat_function(fun=function(x) {(dnorm(x, mean=mean, sd=sd)) * n * binwidth * lambda}, color="black")}, 
  mean=normalmix[[1]]$mu, sd=normalmix[[1]]$sigma, lambda=normalmix[[1]]$lambda, n=length(peak.subclass.TEfrac$ratio),
  binwidth=(max(peak.subclass.TEfrac$ratio)-min(peak.subclass.TEfrac$ratio))/30 #binwidth used for histogram
)
cutoff <- normalmix[["cutoff1"]]
p.val <- pnorm(q=cutoff, mean=normalmix[[1]]$mu[1], sd=normalmix[[1]]$sigma[1], lower.tail=FALSE)
pt_hist_fit_anno <- pt_hist_fit + annotate("segment", x=normalmix[["cutoff1"]], xend=normalmix[["cutoff1"]], y=0, yend=15, colour="black", linetype="dashed") + annotate("text", x=normalmix[["cutoff1"]]+0.005, y=12, label="p<0.005") +
  annotate("rect", xmin=0.137, xmax=0.168, ymin=0, ymax=9.5, alpha=.05, color="red", , linetype="dashed") + annotate("text", x=0.15, y=9, label="highTE GLUT")
ggsave(pt_hist_fit_anno, file="./plot/fig5c.TEfrac.cCREs.histogram.stack.fitmodel.pdf", width=12, height=9)
ggsave(pt_hist_TEfrac_facets, file="./plot/fig14b.TEfrac.cCREs.histogram.facets.pdf", width=12, height=9)

# * compare high TEs percentage GLUT
wilcox_groups <- function(TEfrac, test_groups=c("highTE-GLUT", "other-GLUT"), alt="greater"){
  TEfrac %>% filter(groups %in% test_groups) %>% mutate(groups=factor(groups, levels=test_groups)) -> TEfrac.groups 
  test_groups <- wilcox.test(ratio ~ groups, TEfrac.groups, alternative=alt)
  return(test_groups)
}
wilcox_oGLUT = wilcox_groups(peak.subclass.TEfrac, c("highTE-GLUT", "other-GLUT"), alt="greater")
wilcox_GABA = wilcox_groups(peak.subclass.TEfrac, c("highTE-GLUT", "GABA"), alt="greater")
wilcox_NN = wilcox_groups(peak.subclass.TEfrac, c("highTE-GLUT", "NN"), alt="greater")
wilcox_oGLUT_GABA = wilcox_groups(peak.subclass.TEfrac, c("GABA", "other-GLUT"), alt="two.sided")
wilcox_oGLUT_NN = wilcox_groups(peak.subclass.TEfrac, c("NN", "other-GLUT"), alt="two.sided")
wilcox_NN_GABA = wilcox_groups(peak.subclass.TEfrac, c("GABA", "NN"), alt="two.sided")

# TE fraction groups comparison box plot
pt_TEfrac_box <- ggplot(peak.subclass.TEfrac, aes(x=groups, y=ratio, fill=type, alpha=type)) + geom_boxplot() + scale_fill_manual(values=c("#69b3a2", "grey")) + scale_alpha_manual(values=c(1,0.1)) + theme_classic() + 
  theme(legend.position = "none") + xlab("") + ylab("Fraction of open cCREs overlap with TEs in subclass")
pt_TEfrac_boxP <- pt_TEfrac_box + annotate(geom="segment", x=1, xend=2, y=0.17, yend=0.17, colour="black") + annotate(geom="text", x=1.5, y=0.175, label=paste0("p=", signif(wilcox_oGLUT$p.value, 3))) + 
  annotate(geom="segment", x=1, xend=3, y=0.18, yend=0.18, colour="black") + annotate(geom="text", x=2, y=0.185, label=paste0("p=", signif(wilcox_GABA$p.value, 3))) + 
  annotate(geom="segment", x=1, xend=4, y=0.19, yend=0.19, colour="black") + annotate(geom="text", x=2.5, y=0.195, label=paste0("p=", signif(wilcox_NN$p.value, 3))) 
ggsave(pt_TEfrac_boxP, file="./plot/figS14c.TEfrac_groups_boxplot.pdf", width=8, height=8)


# * plot cv.mean
plot_logcpm_TEsumcount_cv.mean <- smScatter_plot(TEs.count.TE.meta.sum, "cv.logcpm", "mean.logcpm", "TE variation across cell subclass", 
    "Coefficient of variation", expression(bold("Avg. log"["2"])*bold(" (CPM+1)")))   

TEs_scatter <- plot_logcpm_TEsumcount_cv.mean + geom_hline(yintercept=0.6, linetype="dashed", color="white") + 
    geom_vline(xintercept=0.36, linetype="dashed", color="white") + annotate(geom="curve", x=0.5, y=10, xend=0.25, yend=9, 
    curvature=.3, arrow=arrow(length=unit(2, "mm")), color="white") + 
    annotate(geom="text", x=0.52, y=10, label=paste0("Invariable TEs\nn = ", nrow(TEs.count.invariableTE.meta)), hjust="left", color="white") +
    annotate(geom="text", x=1, y=4, label=paste0("Variable TEs\nn = ", nrow(TEs.count.variableTE.meta)), hjust="left", color="white")
ggsave("./plot/figS15.TEs_variability.pdf", TEs_scatter, width=8, height=8)




