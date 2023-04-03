#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ComplexHeatmap"))

packdir <- "./package/R"
import::from(.from="colors.R", .directory=packdir, ArchRPalettes)

# * meta functions
cal_cv <- function(x){
  cv <- sd(x) / mean(x)
  return(cv)
}

meta.TE.sumcount <- function(anno, count, level){
	if(identical(anno$cCRE, colnames(count))){
		anno %>% select(cCRE, AnnoDetail, AnnoFamily, AnnoSuperFam) %>% cbind(t(count)) -> TEs.meta
		colnames(TEs.meta) <- gsub("-", ".", colnames(TEs.meta))
		TEs.meta %>% gather(Subclass, count, rownames(count)) -> TEs.meta.long
		TEs.meta.long %>% group_by_(level, "Subclass") %>% summarise(num=n(), counts=sum(count)) %>% 
			mutate(CPM=counts/total.count.subclass*1000000, logCPM=log2(CPM+1)) -> TEs.meta.sum
		TEs.meta.sum %>% group_by_(level) %>% summarise(nums=n(), mean.cpm=mean(CPM), cv.cpm=cal_cv(CPM), 
			median.cpm=median(CPM), max.cpm=max(CPM), min.cpm=min(CPM), delta.cpm=max(CPM)-min(CPM),
			mean.logcpm=mean(logCPM), cv.logcpm=cal_cv(logCPM), median.logcpm=median(logCPM), 
			max.logcpm=max(logCPM), min.logcpm=min(logCPM), delta.logcpm=max(logCPM)-min(logCPM)) -> TEs.meta.sum.stats
		return(TEs.meta.sum.stats)
	}
}

peak.num <- function(anno, TE.meta){
	anno %>% filter(AnnoDetail %in% TE.meta$AnnoDetail) -> anno.TE.select
	anno.TE.select %>% select(cCRE) %>% distinct() %>% summarise(num=n()) -> TE.peak.num
	return(TE.peak.num)
}

peak.num.set <- function(anno, TE.meta){
	anno %>% filter(AnnoDetail %in% TE.meta$AnnoDetail) -> anno.TE.select
	anno.TE.select %>% group_by(AnnoDetail) %>% summarise(num=n()) -> TEs.peak.num.set
	return(TEs.peak.num.set)
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
peak.anno <- readRDS("./CEMBA_data/rds/update_pmat/peak.anno.rds")
peak.anno.TEs <- readRDS("./CEMBA_data/rds/update_pmat/peak.anno.TEs.rds")
peak.count <- readRDS("./CEMBA_data/rds/update_pmat/mba.whole.subclass.pmat.sc.filter.totalcount.rds")
if(identical(peak.anno$Position, colnames(peak.count))){
	colnames(peak.count) <- peak.anno$cCRE
}
total.count.subclass <- apply(peak.count, 1, sum)
total.count.subclass <- total.count.subclass[order(names(total.count.subclass))]


# * TEs CPM and variablity
TEs.count.TE.meta.sum <- meta.TE.sumcount(peak.anno.TEs, peak.count.TEs, "AnnoDetail")
plot_logcpm_TEsumcount_cv.mean <- smScatter_plot(TEs.count.TE.meta.sum, "cv.logcpm", "mean.logcpm", "TE variation across cell subclass", 
	"Coefficient of variation", expression(bold("Avg. log"["2"])*bold(" (CPM+1)")))   

TEs.count.TE.meta.sum %>% filter(cv.logcpm<0.36, mean.logcpm>0.6) -> TEs.count.invariableTE.meta    
TEs.count.TE.meta.sum %>% filter(cv.logcpm>=0.36, mean.logcpm>0.6) -> TEs.count.variableTE.meta 
TEvariable.peak.num <- peak.num(peak.anno.TEs, TEs.count.variableTE.meta)
TEinvariable.peak.num <- peak.num(peak.anno.TEs, TEs.count.invariableTE.meta)
TEvariable.peak.num.set <- peak.num.set(peak.anno.TEs, TEs.count.variableTE.meta)
TEinvariable.peak.num.set <- peak.num.set(peak.anno.TEs, TEs.count.invariableTE.meta)

# * plot output
TEs_scatter <- plot_logcpm_TEsumcount_cv.mean + geom_hline(yintercept=0.6, linetype="dashed", color="white") + 
	geom_vline(xintercept=0.36, linetype="dashed", color="white") + annotate(geom="curve", x=0.5, y=10, xend=0.25, yend=9, 
    curvature=.3, arrow=arrow(length=unit(2, "mm")), color="white") + 
    annotate(geom="text", x=0.52, y=10, label=paste0("Invariable TEs\nn = ", nrow(TEs.count.invariableTE.meta)), hjust="left", color="white") +
    annotate(geom="text", x=1, y=4, label=paste0("Variable TEs\nn = ", nrow(TEs.count.variableTE.meta)), hjust="left", color="white")
ggsave("./plot/figS15.TEs_variability.pdf", TEs_scatter, width=8, height=8)


