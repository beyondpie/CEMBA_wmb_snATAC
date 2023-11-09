suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))

# * functions from Songpeng's package
projdir <- here::here()
packdir <- file.path(projdir, "scatac/scatac_mousebrain/package/R")
import::from(.from = "ggtheme.R", .directory = packdir, dotplotTheme)
import::from(.from="colors.R", .directory=packdir, SnapATACPalette)

getMaxColScoreWithName <- function(mat) {
  maxScores <- apply(mat, 2, max)
  rnames <- rownames(mat)
  maxNames <- apply(mat, 2, \(x){
    rnames[which.max(x)]
  })
  names(maxScores) <- maxNames
  return(maxScores)
}
to3.matrix <- function(mat,
                       names = c("row", "column", "score"),
                       int2str = TRUE,
                       factor2str = TRUE) {
  r <- reshape2::melt(mat)
  colnames(r) <- names
  if (factor2str) {
    if (is.factor(r[,1])) {
      r[,1] <- as.character(r[,1])
    }
    if (is.factor(r[,2])) {
      r[,2] <- as.character(r[,2])
    }
  }
  if(int2str) {
    if(is.numeric(r[,1])) {
      r[,1] <- as.character(r[,1])
    }
    if(is.numeric(r[,2])) {
      r[,2] <- as.character(r[,2])
    }
  }
  return(r)
}
prepareDotPlot4TransferLabel <- function(tfmat,
                                         refOrder = NULL,
                                         names = c("row", "column", "score"),
                                         ignoreEmptyRef = TRUE) {
  maxScore <- getMaxColScoreWithName(mat = tfmat)
  query2ref <- data.frame(
    query = colnames(tfmat),
    ref = names(maxScore),
    row.names = colnames(tfmat)
  )
  if (ignoreEmptyRef) {
    message("remove refs not having query mapped to.")
    tfmat <- tfmat[rownames(tfmat) %in% query2ref$ref, ]
  }
  if (is.null(refOrder)) {
    message("refOrder is null, will use default numeric order for it.")
    refOrder <- rownames(tfmat)
    refOrder <- refOrder[order(as.integer(refOrder))]
  } else {
    refOrder <- refOrder[refOrder %in% rownames(tfmat)]
  }
  queryOrder <- query2ref$query[
    order(factor(query2ref$ref, levels = refOrder))]

  meltmat <- to3.matrix(tfmat, names)
  meltmat[,1] <- factor(meltmat[,1], levels = refOrder)
  meltmat[,2] <- factor(meltmat[,2], levels = queryOrder)
  # reduce size of meltmat
  meltmat <- meltmat[meltmat[,3] > 0, ]
  return(meltmat)
}

# * load data
meta_v9.7 <- readRDS("mba.whole.cell.meta.v9.7.rds")
all_score <- read.csv("sa2.all.subclass2L4.transferLabelScore.csv")
allen_dict <- read.csv("allen_subclass.dict.csv")
atac_dict <- read.csv("atac.subclass.dict.csv")

# * split neuron and non-neuron
meta_v9.7 %>% select(L4, subclass_id_v3, NT_v3) %>% distinct() -> L4_subclass.dict
L4_subclass.dict %>% filter(NT_v3 == "NN") %>% mutate(L4=paste0("X", L4)) %>% mutate(L4=gsub("-",".",L4)) -> L4_subclass.dict_nn
L4_subclass.dict %>% filter(NT_v3 != "NN") %>% mutate(L4=paste0("X", L4)) %>% mutate(L4=gsub("-",".",L4)) -> L4_subclass.dict_neu
meta_v9.7_neu <- meta_v9.7 %>% dplyr::filter(NT_v3 != "NN")
meta_v9.7_nn <- meta_v9.7 %>% dplyr::filter(NT_v3 == "NN")

# * 253 neuronal subclass score
all_score <- all_score[match(allen_dict$subclass_id, rownames(all_score)),]
#rownames(all_score) <- paste0("allen.", rownames(all_score))
all_score <- as.matrix(all_score)
neu253_score <- all_score[match(atac_dict$subclass_id[which(atac_dict$NT != "NN")], rownames(all_score)), match(L4_subclass.dict_neu$L4, colnames(all_score))]

mat.dotplot_neu253 <- prepareDotPlot4TransferLabel(tfmat = neu253_score, ignoreEmptyRef = FALSE, names = c("Allen", "ATAC", "score"))
allen.miss <- setdiff(rownames(neu253_score), mat.dotplot_neu253$Allen)
#mat.dotplot_neu253.miss <- data.frame(Allen = allen.miss, ATAC = mat.dotplot_neu253$ATAC[1], score = 0)
#mat.dotplot_neu253.all <- rbind(mat.dotplot_neu253, mat.dotplot_neu253.miss)

mytheme <- dotplotTheme(legend.pos = "right")

lowSimScore<- quantile(mat.dotplot_neu253$score, 0.1)
#highSimScore <- max(mat.dotplot_neu253$score)
highSimScore <- quantile(mat.dotplot_neu253$score, 1)

palette <- "Reds"
p.tfvote <- ggplot(data = mat.dotplot_neu253, aes(x = ATAC, y = Allen)) +
  geom_tile(aes(fill = score)) +
  mytheme +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()) +
  #scale_fill_gradient(low = "white", high = "red",
  #  limits = c(lowSimScore, highSimScore), na.value = "white") +
  scale_fill_distiller(palette = palette, direction = 1, limits = c(lowSimScore, highSimScore), na.value = "white") + 
  #scale_size(range = c(0, 1)) +
  xlab(paste0(ncol(neu253_score), " clusters of snATAC-seq")) +
  ylab(paste0(nrow(neu253_score), " clusters of scRNA-seq")) 
ggsave(p.tfvote, file="Fig1e.subclass.neu253.consensus_score.tile.red.pdf", width=16, height=16)

# * plot1 class
meta_v9.7_neu %>% select(subclass_id_label_v3, class_id_label_v3, class_color_v3) %>% distinct() %>% arrange(subclass_id_label_v3) %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)), 
		class_id_label_v3=factor(class_id_label_v3, levels=unique(class_id_label_v3)),
		class_color_v3=factor(class_color_v3, levels=unique(class_color_v3))) -> meta_v9.7_neu_subclass.class
pt_tile_neu_class <- ggplot(meta_v9.7_neu_subclass.class, aes(x=0.2, y=subclass_id_label_v3, fill=class_id_label_v3)) + 
	geom_tile() + 
	scale_fill_manual(values=levels(meta_v9.7_neu_subclass.class$class_color_v3)) + scale_x_discrete(expand=expansion(0)) +  
	theme_minimal() + xlab("Class") + ylab("253 non-neuronal subclass") + guides(fill=guide_legend(title="Class")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_blank(), axis.text.y=element_text(family="serif", size=4, color="black"), axis.title.y=element_text(family="serif", color="black"), 
		axis.title.x=element_text(family="serif", angle=90, vjust=0.5, hjust=1, color="black")) 
meta_v9.7_nn %>% select(subclass_id_label_v3, class_id_label_v3, class_color_v3) %>% distinct() %>% arrange(subclass_id_label_v3) %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)), 
		class_id_label_v3=factor(class_id_label_v3, levels=unique(class_id_label_v3)),
		class_color_v3=factor(class_color_v3, levels=unique(class_color_v3))) -> meta_v9.7_nn_subclass.class
pt_tile_nn_class <- ggplot(meta_v9.7_nn_subclass.class, aes(x=0.2, y=subclass_id_label_v3, fill=class_id_label_v3)) + 
	geom_tile() + 
	scale_fill_manual(values=levels(meta_v9.7_nn_subclass.class$class_color_v3)) + scale_x_discrete(expand=expansion(0)) +  
	theme_minimal() + xlab("Class") + ylab("22 non-neuronal subclass") + guides(fill=guide_legend(title="Class")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_blank(), axis.text.y=element_text(family="serif", size=4, color="black"), axis.title.y=element_text(family="serif", color="black"), 
		axis.title.x=element_text(family="serif", angle=90, vjust=0.5, hjust=1, color="black")) 
# * plot2 NT
meta_v9.7_neu %>% select(subclass_id_label_v3, NT_v3, nt_type_color_v3) %>% distinct(subclass_id_label_v3, .keep_all=TRUE) %>% arrange(subclass_id_label_v3) %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)), 
		NT_v3=factor(NT_v3, levels=unique(NT_v3)), 
		nt_type_color_v3=factor(nt_type_color_v3, levels=unique(nt_type_color_v3)))-> meta_v9.7_neu_subclass.nt
pt_tile_neu_nt <- ggplot(meta_v9.7_neu_subclass.nt, aes(x=0.2, y=subclass_id_label_v3, fill=NT_v3)) + 
	geom_tile() + 
	scale_fill_manual(values=levels(meta_v9.7_neu_subclass.nt$nt_type_color_v3)) + scale_x_discrete(expand=expansion(0)) +  
	theme_minimal() + xlab("NT type") + guides(fill=guide_legend(title="NT type")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), 
		axis.title.x=element_text(family="serif", angle=90, vjust=0.5, hjust=1, color="black"))
# * plot3 replicate
meta_v9.7_neu %>% select(subclass_id_label_v3, biorep) %>% arrange(subclass_id_label_v3) %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3))) -> meta_v9.7_neu_subclass_biorep
pt_bar_neu_biorep <- ggplot(meta_v9.7_neu_subclass_biorep, aes(subclass_id_label_v3, fill = biorep)) + geom_bar(position = "fill") + 
	scale_fill_manual(labels=c("rep1", "rep2"), values=brewer.pal(length(unique(meta_v9.7_neu_subclass_biorep$biorep)), "Paired")) + scale_y_discrete(expand=expansion(0)) + 
	theme_minimal() + ylab("Replicate") + xlab("") + guides(fill=guide_legend(title="Replicate")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_text(family="serif", color="black"), axis.title.x=element_text(family="serif", color="black")) + 
	coord_flip()
meta_v9.7_nn %>% select(subclass_id_label_v3, biorep) %>% arrange(subclass_id_label_v3) %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3))) -> meta_v9.7_nn_subclass_biorep
pt_bar_nn_biorep <- ggplot(meta_v9.7_nn_subclass_biorep, aes(subclass_id_label_v3, fill = biorep)) + geom_bar(position = "fill") + 
	scale_fill_manual(labels=c("rep1", "rep2"), values=brewer.pal(length(unique(meta_v9.7_nn_subclass_biorep$biorep)), "Paired")) + scale_y_discrete(expand=expansion(0)) + 
	theme_minimal() + ylab("Replicate") + xlab("") + guides(fill=guide_legend(title="Replicate")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_text(family="serif", color="black"), axis.title.x=element_text(family="serif", color="black")) + 
	coord_flip()
# * plot4 region
meta_v9.7_neu %>% select(MajorRegion, MajorRegionColor) %>% distinct() %>% arrange(MajorRegion) -> region.dict
meta_v9.7_neu %>% select(subclass_id_label_v3, MajorRegion, MajorRegionColor) %>% arrange(subclass_id_label_v3) %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)),
		MajorRegion=factor(MajorRegion, levels=region.dict$MajorRegion),
		MajorRegionColor=factor(MajorRegionColor, levels=region.dict$MajorRegionColor)) -> meta_v9.7_neu_subclass_region
pt_bar_neu_region <- ggplot(meta_v9.7_neu_subclass_region, aes(subclass_id_label_v3, fill = MajorRegion)) + geom_bar(position = "fill") + 
	scale_fill_manual(values=levels(meta_v9.7_neu_subclass_region$MajorRegionColor), labels=gsub("Isocortex", "CTX", levels(meta_v9.7_neu_subclass_region$MajorRegion))) + 
	scale_y_discrete(expand=expansion(0)) + 
	theme_minimal() + ylab("Major Region") + xlab("") + guides(fill=guide_legend(title="Major Region")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_text(family="serif", color="black"), 
		axis.title.x=element_text(family="serif", color="black")) + 
	coord_flip()
meta_v9.7_nn %>% select(subclass_id_label_v3, MajorRegion, MajorRegionColor) %>% arrange(subclass_id_label_v3) %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)),
		MajorRegion=factor(MajorRegion, levels=region.dict$MajorRegion),
		MajorRegionColor=factor(MajorRegionColor, levels=region.dict$MajorRegionColor)) -> meta_v9.7_nn_subclass_region
pt_bar_nn_region <- ggplot(meta_v9.7_nn_subclass_region, aes(subclass_id_label_v3, fill = MajorRegion)) + geom_bar(position = "fill") + 
	scale_fill_manual(values=levels(meta_v9.7_nn_subclass_region$MajorRegionColor), labels=gsub("Isocortex", "CTX", levels(meta_v9.7_nn_subclass_region$MajorRegion))) + 
	scale_y_discrete(expand=expansion(0)) + 
	theme_minimal() + ylab("Major Region") + xlab("") + guides(fill=guide_legend(title="Major Region")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_text(family="serif", color="black"), 
		axis.title.x=element_text(family="serif", color="black")) + 
	coord_flip()
# * plot5 ncluster
meta_v9.7_neu %>% select(subclass_id_label_v3, class_id_label_v3, class_color_v3, L4) %>% distinct() %>% 
	group_by(subclass_id_label_v3, class_id_label_v3, class_color_v3) %>% summarise(nL4=n()) %>% ungroup %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)), 
		class_id_label_v3=factor(class_id_label_v3, levels=unique(class_id_label_v3)),
		class_color_v3=factor(class_color_v3, levels=unique(class_color_v3))) -> meta_v9.7_neu_subclass.L4s
#breaks <- 2^c(log2(1), log2(5), log2(10), log2(20), log2(30), log2(40)) 
pt_bar_neu_L4 <- ggplot(meta_v9.7_neu_subclass.L4s, aes(x=subclass_id_label_v3, y=nL4, fill = class_id_label_v3)) + geom_col() + 
	scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=2), breaks=c(1,5,10,20,30,40), expand=expansion(0)) + 
	geom_hline(yintercept=1, linetype="dashed", color="grey") + geom_hline(yintercept=5, linetype="dashed", color="grey") + geom_hline(yintercept=10, linetype="dashed", color="grey") + 
	geom_hline(yintercept=20, linetype="dashed", color="grey") +geom_hline(yintercept=30, linetype="dashed", color="grey") + geom_hline(yintercept=40, linetype="dashed", color="grey") + 
	scale_fill_manual(values=levels(meta_v9.7_neu_subclass.L4s$class_color_v3)) + 
	theme_minimal() + ylab("# of clusters per subclass") + xlab("") + guides(fill=guide_legend(title="Class")) +
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_text(family="serif", angle=90, vjust=0.5, hjust=1, color="black"), axis.text.y=element_blank(), axis.title.y=element_text(family="serif", color="black"), axis.title.x=element_text(family="serif", color="black")) + 
	coord_flip()
meta_v9.7_nn %>% select(subclass_id_label_v3, class_id_label_v3, class_color_v3, L4) %>% distinct() %>% 
	group_by(subclass_id_label_v3, class_id_label_v3, class_color_v3) %>% summarise(nL4=n()) %>% ungroup %>% 
	mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)), 
		class_id_label_v3=factor(class_id_label_v3, levels=unique(class_id_label_v3)),
		class_color_v3=factor(class_color_v3, levels=unique(class_color_v3))) -> meta_v9.7_nn_subclass.L4s
#breaks <- 2^c(log2(1), log2(5), log2(10), log2(20), log2(30), log2(40)) 
pt_bar_nn_L4 <- ggplot(meta_v9.7_nn_subclass.L4s, aes(x=subclass_id_label_v3, y=nL4, fill = class_id_label_v3)) + geom_col() + 
	scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=2), breaks=c(1,5,10,20,30,40,50), expand=expansion(0)) + 
	geom_hline(yintercept=1, linetype="dashed", color="grey") + geom_hline(yintercept=5, linetype="dashed", color="grey") + geom_hline(yintercept=10, linetype="dashed", color="grey") + 
	geom_hline(yintercept=20, linetype="dashed", color="grey") +geom_hline(yintercept=30, linetype="dashed", color="grey") + geom_hline(yintercept=40, linetype="dashed", color="grey") + 
	scale_fill_manual(values=levels(meta_v9.7_nn_subclass.L4s$class_color_v3)) + 
	theme_minimal() + ylab("# of clusters per subclass") + xlab("") + guides(fill=guide_legend(title="Class")) +
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_text(family="serif", angle=90, vjust=0.5, hjust=1, color="black"), axis.text.y=element_blank(), axis.title.y=element_text(family="serif", color="black"), axis.title.x=element_text(family="serif", color="black")) + 
	coord_flip()
# * plot6 ncell
meta_v9.7_neu %>% group_by(subclass_id_label_v3, class_id_label_v3, class_color_v3) %>% summarise(ncell=n()) %>% 
	ungroup %>% mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)), 
	class_id_label_v3=factor(class_id_label_v3, levels=unique(class_id_label_v3)),
	class_color_v3=factor(class_color_v3, levels=unique(class_color_v3))) -> meta_v9.7_neu_subclass.ncells
pt_bar_neu_ncell <- ggplot(meta_v9.7_neu_subclass.ncells, aes(x=subclass_id_label_v3, y=ncell, fill = class_id_label_v3)) + geom_col() + 
	scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=10), breaks=(10**c(1:5)), labels=scales::trans_format('log10',scales::math_format(10^.x)), expand=expansion(0)) + 
	geom_hline(yintercept=10, linetype="dashed", color="grey") + geom_hline(yintercept=100, linetype="dashed", color="grey") + geom_hline(yintercept=1000, linetype="dashed", color="grey") + 
	geom_hline(yintercept=10000, linetype="dashed", color="grey") + geom_hline(yintercept=100000, linetype="dashed", color="grey") + 
	scale_fill_manual(values=levels(meta_v9.7_neu_subclass.ncells$class_color_v3)) +  
	theme_minimal() + ylab("# of nuclei per subclass") + xlab("") + guides(fill=guide_legend(title="Class")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_text(family="serif", angle=90, vjust=0.5, hjust=1, color="black"), axis.text.y=element_blank(), axis.title.y=element_text(family="serif", color="black"), axis.title.x=element_text(family="serif", color="black")) + 
	coord_flip()
meta_v9.7_nn %>% group_by(subclass_id_label_v3, class_id_label_v3, class_color_v3) %>% summarise(ncell=n()) %>% 
	ungroup %>% mutate(subclass_id_label_v3=factor(subclass_id_label_v3, levels=unique(subclass_id_label_v3)), 
	class_id_label_v3=factor(class_id_label_v3, levels=unique(class_id_label_v3)),
	class_color_v3=factor(class_color_v3, levels=unique(class_color_v3))) -> meta_v9.7_nn_subclass.ncells
pt_bar_nn_ncell <- ggplot(meta_v9.7_nn_subclass.ncells, aes(x=subclass_id_label_v3, y=ncell, fill = class_id_label_v3)) + geom_col() + 
	scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1, base=10), breaks=(10**c(1:5)), labels=scales::trans_format('log10',scales::math_format(10^.x)), expand=expansion(0)) + 
	geom_hline(yintercept=10, linetype="dashed", color="grey") + geom_hline(yintercept=100, linetype="dashed", color="grey") + geom_hline(yintercept=1000, linetype="dashed", color="grey") + 
	geom_hline(yintercept=10000, linetype="dashed", color="grey") + geom_hline(yintercept=100000, linetype="dashed", color="grey") + 
	scale_fill_manual(values=levels(meta_v9.7_nn_subclass.ncells$class_color_v3)) +  
	theme_minimal() + ylab("# of nuclei per subclass") + xlab("") + guides(fill=guide_legend(title="Class")) + 
	theme(legend.position="bottom", legend.title=element_text(family="serif", size=8, angle=90, vjust=0.5, hjust=1), legend.text=element_text(family="serif", size=6), 
		axis.text.x=element_text(family="serif", angle=90, vjust=0.5, hjust=1, color="black"), axis.text.y=element_blank(), axis.title.y=element_text(family="serif", color="black"), axis.title.x=element_text(family="serif", color="black")) + 
	coord_flip()
# * combined plots
pt_neu_subclass_bar <- ggarrange(pt_tile_neu_class, pt_tile_neu_nt, pt_bar_neu_biorep, pt_bar_neu_region, pt_bar_neu_L4, pt_bar_neu_ncell,
          nrow=1,  align="h", 
          widths=c(10, 5, 20, 20, 20, 20))
pt_nn_subclass_bar <- ggarrange(pt_tile_nn_class, pt_bar_nn_biorep, pt_bar_nn_region, pt_bar_nn_L4, pt_bar_nn_ncell,
          nrow=1,  align="h", 
          widths=c(10, 20, 20, 20, 20))
ggsave(file="./fig1f_pt_neu_subclass_bar.pdf", pt_neu_subclass_bar, width=20, height=32)
ggsave(file="./figS6i_pt_nn_subclass_bar.pdf", pt_nn_subclass_bar, width=20, height=15)

