#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("mixtools"))

# * load data 
peak.anno <- readRDS("./CEMBA_data/rds/update_pmat/peak.anno.rds")
peak.subclass <- readRDS("./CEMBA_data/rds/update_pmat/peak.subclass.rds")
peak.anno.TEs <- readRDS("./CEMBA_data/rds/update_pmat/peak.anno.TEs.rds")
celltype.dict <- readRDS("./CEMBA_data/rds/update_pmat/celltype.dict.rds")
L2.order <- readRDS("./CEMBA_data/rds/update_pmat/L2.order.rds")
peak.subclass.TEs <- peak.subclass %>% as.data.frame %>% select(peak.anno.TEs$cCRE)
identical(colnames(peak.subclass.TEs), peak.anno.TEs$cCRE) #TRUE

# * functions
concatAllenAnnot1 <- function(allenAnnot, space2="_", slash2="-") {
  r <- gsub("/", slash2, allenAnnot)
  r <- gsub(" ", space2, r)
  r <- gsub("-", ".", r)
  return(r)}

meta.peak.class <- function(anno, pmat, level){
	if(identical(colnames(pmat),anno$cCRE)){
		rownames(pmat) <- concatAllenAnnot1(rownames(pmat))
		anno %>% select(cCRE, AnnoDetail, AnnoFamily1, AnnoSuperFam1) %>% cbind(t(pmat)) -> peak.meta
		peak.meta %>% gather(Subclass, peakmat, rownames(pmat)) -> peak.meta.long
		peak.meta.long %>% group_by_("Subclass", level) %>% summarise(num=n(), counts=sum(peakmat)) %>%
			mutate(ratio=counts/sum(counts)) -> peak.meta.subclass
		#mutate(label=paste(level, scales::percent(as.numeric(ratio), 0.1),sep=', ')) %>% mutate(Ratio=str_replace(label,", 0.0%", ""))
		peak.meta.subclass %>% select(Subclass, level, ratio) %>% spread(Subclass, ratio) %>% column_to_rownames(level) -> peak.meta.subclass.ratio
		return(peak.meta.subclass.ratio)
	}
}

# * subclass cCRE annotation heatmap
peak.subclass.AnnoRatio <- meta.peak.class(peak.anno, peak.subclass, "AnnoSuperFam1")
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
peak.subclass %>% t %>% as.data.frame %>% rownames_to_column(var="cCRE") -> peak.subclass.data
peak.subclass.data %>% mutate(TEs=ifelse(cCRE %in% peak.anno.TEs$cCRE,"TEs","nonTEs")) %>% 
	gather(Subclass, peakmat, rownames(peak.subclass)) -> peak.subclass.long 
peak.subclass.long %>% group_by(Subclass, TEs) %>% summarise(num=n(), counts=sum(peakmat)) %>% 
	mutate(ratio=counts/sum(counts)) %>% filter(TEs=="TEs")	-> peak.subclass.TEfrac
peak.subclass.TEfrac$mainclass <- celltype.dict$mainclass[match(peak.subclass.TEfrac$Subclass, celltype.dict$AllenAnnotConcat)]
#peak.subclass.TEfrac$class <- celltype.dict$class[match(peak.subclass.TEfrac$Subclass, celltype.dict$AllenAnnotConcat)]

cols_mainclass <- celltype.dict %>% select(mainclass, mainclassColorv2) %>% distinct %>% select(mainclassColorv2) %>% unlist
names(cols_mainclass) <- celltype.dict %>% select(mainclass, mainclassColorv2) %>% distinct %>% select(mainclass) %>% unlist

pt_hist_TEfrac <- ggplot(peak.subclass.TEfrac, aes(x=ratio, fill=mainclass, color=mainclass)) + geom_histogram(alpha=0.6, position="identity") +
	theme_classic() + scale_color_manual(values=cols_mainclass) + scale_fill_manual(values=cols_mainclass) + labs(x="Fraction of cCREs overlpped with TEs in subclass", y="# of cell subclasses") + 
	annotate("text", x=0.05, y=11, label="Cell subclass")

pt_hist_TEfrac_facets <- ggplot(peak.subclass.TEfrac, aes(x=ratio)) + geom_histogram(aes(color=mainclass, fill=mainclass)) + facet_grid(mainclass ~ .) +
	theme_classic() + scale_color_manual(values=cols_mainclass) + scale_fill_manual(values=cols_mainclass) + theme(legend.position="none") + 
	labs(x="Fraction of cCREs overlpped with TEs in subclass", y="# of cell subclasses")

# * highTE cell subclasses (mixture model)
fit_model_auto <- function(TEfrac, epsilon=1e-03, p1=0.005, p2=0.5, maxit=10000){
  set.seed(2023)
  normalmix <- normalmixEM(TEfrac$ratio,k=2,epsilon=epsilon, maxit=maxit)
  pi <- normalmix$lambda[1]
  mu1 <- normalmix$mu[1]
  mu2 <- normalmix$mu[2]
  sigma1 <- normalmix$sigma[1]
  sigma2 <- normalmix$sigma[2]
#  y1 <- dnorm(x, mu1, sigma1)
#  y2 <- dnorm(x, mu2, sigma2)
  Px1 <- qnorm(p=p1, mean=mu1, sd=sigma1, lower.tail=FALSE, log.p = FALSE)
  Px2 <- qnorm(p=p2, mean=mu2, sd=sigma2, lower.tail=FALSE, log.p = FALSE)
  cutoff1 <- round(Px1, 3)
  cutoff2 <- round(Px2, 3)
  NormIntersect <- function(x, mu1, sig1, mu2, sig2) {
    (dnorm(x, mean=mu1, sd=sig1) - dnorm(x, mean=mu2, sd=sig2))**2
  }
  Intersect <- optimize(NormIntersect, interval=c(mu1,mu2), mu1=mu1, sig1=sigma1, 
                      mu2=mu2, sig2=sigma2)
  intersectP <- Intersect$minimum
  cutoff3 <- round(intersectP,3)
  return(fit_EM=list(normalmix=normalmix, pi=pi, mu1=mu1, mu2=mu2, sigma1=sigma1, sigma2=sigma2, cutoff1=cutoff1, cutoff2=cutoff2, cutoff3=cutoff3))
} 

normalmix <- fit_model_auto(peak.subclass.TEfrac, epsilon=1e-18, p1=0.005, p2=0.5)
pt_hist_TEfrac <- ggplot(peak.subclass.TEfrac, aes(x=ratio, fill=mainclass, color=mainclass)) + geom_histogram(alpha=0.6, position="stack", binwidth=(max(peak.subclass.TEfrac$ratio)-min(peak.subclass.TEfrac$ratio))/150) +
  theme_classic() +scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2") + labs(x="Fraction of open cCREs overlap with TEs in subclass", y="Number of cell subclasses") + 
  annotate("text", x=0.05, y=11, label="Cell subclass")
pt_hist_fit <- pt_hist_TEfrac + mapply(function(mean, sd, lambda, n, binwidth) {
  stat_function(fun=function(x) {(dnorm(x, mean=mean, sd=sd)) * n * binwidth * lambda}, color="black")}, 
  mean=normalmix[[1]]$mu, sd=normalmix[[1]]$sigma, lambda=normalmix[[1]]$lambda, n=length(peak.subclass.TEfrac$ratio),
  binwidth=(max(peak.subclass.TEfrac$ratio)-min(peak.subclass.TEfrac$ratio))/30 #binwidth used for histogram
)
p.val <- pnorm(q=0.137, mean=normalmix[[1]]$mu[1], sd=normalmix[[1]]$sigma[1], lower.tail=FALSE)
pt_hist_fit_anno <- pt_hist_fit + annotate("segment", x=normalmix[["cutoff1"]], xend=normalmix[["cutoff1"]], y=0, yend=15, colour="black", linetype="dashed") + annotate("text", x=normalmix[["cutoff1"]]+0.005, y=12, label="p<0.005") +
  annotate("rect", xmin=0.137, xmax=0.168, ymin=0, ymax=9.5, alpha=.05, color="red", , linetype="dashed") + annotate("text", x=0.15, y=9, label="highTE GLUT")
ggsave(pt_hist_fit_anno, file="./plot/fig5c.TEfrac.cCREs.histogram.stack.fitmodel.pdf", width=12, height=9)
ggsave(pt_hist_TEfrac_facets, file="./plot/fig14b.TEfrac.cCREs.histogram.facets.pdf", width=12, height=9)

# * compare high TEs percentage GLUT
peak.subclass.TEfrac %>% arrange(desc(ratio)) %>% print(n=20)
peak.subclass.TEfrac %>% mutate(groups=if_else(ratio>0.137 & mainclass=="GLUT", "highTE-GLUT", if_else(ratio<=0.137 & mainclass=="GLUT", "other-GLUT", mainclass))) -> peak.subclass.TEfrac
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
peak.subclass.TEfrac %>% mutate(groups=factor(groups, levels=c("highTE-GLUT", "other-GLUT", "GABA", "NN"))) %>% mutate(type=ifelse(groups=="highTE-GLUT","Highlighted","Normal")) -> peak.subclass.TEfrac
pt_TEfrac_box <- ggplot(peak.subclass.TEfrac, aes(x=groups, y=ratio, fill=type, alpha=type)) + geom_boxplot() + scale_fill_manual(values=c("#69b3a2", "grey")) + scale_alpha_manual(values=c(1,0.1)) + theme_classic() + 
  theme(legend.position = "none") + xlab("") + ylab("Fraction of open cCREs overlap with TEs in subclass")
pt_TEfrac_boxP <- pt_TEfrac_box + annotate(geom="segment", x=1, xend=2, y=0.17, yend=0.17, colour="black") + annotate(geom="text", x=1.5, y=0.175, label=paste0("p=", signif(wilcox_oGLUT$p.value, 3))) + 
  annotate(geom="segment", x=1, xend=3, y=0.18, yend=0.18, colour="black") + annotate(geom="text", x=2, y=0.185, label=paste0("p=", signif(wilcox_GABA$p.value, 3))) + 
  annotate(geom="segment", x=1, xend=4, y=0.19, yend=0.19, colour="black") + annotate(geom="text", x=2.5, y=0.195, label=paste0("p=", signif(wilcox_NN$p.value, 3))) 
ggsave(pt_TEfrac_boxP, file="./plot/figS14c.TEfrac_groups_boxplot.pdf", width=8, height=8)
