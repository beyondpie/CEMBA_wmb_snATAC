#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))


# * function
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
peak.subclass.CAcons.sumSubclass %>% as.data.frame %>% rownames_to_column %>% rename(Subclasses=".") -> peak.subclass.CAcons.sumSubclass
peak.subclass.CAdiver.sumSubclass %>% as.data.frame %>% rownames_to_column %>% rename(Subclasses=".") -> peak.subclass.CAdiver.sumSubclass

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

# * save data
saveRDS(Anno_sepcVSortho_totalcount, "./CEMBA_data/rds/Anno_sepcVSortho_totalcount.rds")
saveRDS(peak.specific.AnnoSuper, "./CEMBA_data/rds/peak.specific.AnnoSuper.rds")
saveRDS(Anno_specVSortho_count, "./CEMBA_data/rds/Anno_specVSortho_count.rds")
saveRDS(Anno_consVSdiver_totalcount, "./CEMBA_data/rds/Anno_consVSdiver_totalcount.rds")
saveRDS(Anno_consVSdiver_ratio, "./CEMBA_data/rds/Anno_consVSdiver_ratio.rds")
saveRDS(Anno_specVSortho_ratio, "./CEMBA_data/rds/Anno_specVSortho_ratio.rds")
saveRDS(sumPeak_CAcons_CAdiver, "./CEMBA_data/rds/sumPeak_CAcons_CAdiver.rds")
saveRDS(peak.subclass.CAcons.sumSubclass, "./CEMBA_data/rds/peak.subclass.CAcons.sumSubclass.rds")
saveRDS(peak.subclass.CAdiver.sumSubclass, "./CEMBA_data/rds/peak.subclass.CAdiver.sumSubclass.rds")
saveRDS(peak.subclass.data.long.fraction_cons.diver, "./CEMBA_data/rds/peak.subclass.data.long.fraction_cons.diver.rds")





