#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))


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

TEs.count.TE.meta.sum %>% filter(cv.logcpm<0.36, mean.logcpm>0.6) -> TEs.count.invariableTE.meta    
TEs.count.TE.meta.sum %>% filter(cv.logcpm>=0.36, mean.logcpm>0.6) -> TEs.count.variableTE.meta 
TEvariable.peak.num <- peak.num(peak.anno.TEs, TEs.count.variableTE.meta)
TEinvariable.peak.num <- peak.num(peak.anno.TEs, TEs.count.invariableTE.meta)
TEvariable.peak.num.set <- peak.num.set(peak.anno.TEs, TEs.count.variableTE.meta)
TEinvariable.peak.num.set <- peak.num.set(peak.anno.TEs, TEs.count.invariableTE.meta)

# * save data
saveRDS(TEs.count.TE.meta.sum, "TEs.count.TE.meta.sum.rds")
saveRDS(TEs.count.invariableTE.meta, "TEs.count.invariableTE.meta.rds")



