\#!/usr/bin/R

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

# * TE fraction in each subclass histogram plot
peak.subclass %>% t %>% as.data.frame %>% rownames_to_column(var="cCRE") -> peak.subclass.data
peak.subclass.data %>% mutate(TEs=ifelse(cCRE %in% peak.anno.TEs$cCRE,"TEs","nonTEs")) %>% 
	gather(Subclass, peakmat, rownames(peak.subclass)) -> peak.subclass.long 
peak.subclass.long %>% group_by(Subclass, TEs) %>% summarise(num=n(), counts=sum(peakmat)) %>% 
	mutate(ratio=counts/sum(counts)) %>% filter(TEs=="TEs")	-> peak.subclass.TEfrac
peak.subclass.TEfrac$mainclass <- celltype.dict$mainclass[match(peak.subclass.TEfrac$Subclass, celltype.dict$AllenAnnotConcat)]
#peak.subclass.TEfrac$class <- celltype.dict$class[match(peak.subclass.TEfrac$Subclass, celltype.dict$AllenAnnotConcat)]

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

# * compare high TEs percentage GLUT
peak.subclass.TEfrac %>% arrange(desc(ratio)) %>% print(n=20)
peak.subclass.TEfrac %>% mutate(groups=if_else(ratio>0.137 & mainclass=="GLUT", "highTE-GLUT", if_else(ratio<=0.137 & mainclass=="GLUT", "other-GLUT", mainclass))) -> peak.subclass.TEfrac

# TE fraction groups comparison box plot
peak.subclass.TEfrac %>% mutate(groups=factor(groups, levels=c("highTE-GLUT", "other-GLUT", "GABA", "NN"))) %>% mutate(type=ifelse(groups=="highTE-GLUT","Highlighted","Normal")) -> peak.subclass.TEfrac

# * save data
saveRDS(normalmix, "./CEMBA_data/rds/normalmix.rds")
saveRDS(peak.subclass.AnnoRatio, "./CEMBA_data/rds/peak.subclass.AnnoRatio.rds")
saveRDS(peak.subclass.TEfrac, "./CEMBA_data/rds/peak.subclass.TEfrac.rds")

