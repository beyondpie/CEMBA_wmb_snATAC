This folder includes Mouse Whole Brain taxonomy's cellxgene matrix for scRNAseq data (10Xv2, 10Xv3) and snRNAseq data (Multiome) in h5ad file and the cluster annotation file.

# cellxgene matrix

  - filename : AIT21_10Xv2.h5ad
               AIT21_10Xv3.h5ad
               AIT21_10Xmulti.h5ad

  - anndata$X      : log-normalized count matrix
           $layers : raw count matrix
           $obs    : cluster id 'cl' and sample related metadata
                     'cl' - cluster id to be used as matching key for cluster  annotation file!
                     library_prep        - library
                     gene.counts.0       - number of detected genes 
                     doublet_score       - doublet score
                     roi                 - region
                     umi.counts          - number of detected UMI's
                     method              - 10Xv3 / 10Xv2 / 10Xmulti
                     sex                 - 
                     external_donor_name - 
                     age                 -   
                     medical_conditions  - Light/Dark

# cluster annotation 

  - filename : AIT21_annotation.tsv
        cl                              : cluster key matching all cluster-related tables
	cluster_id/cluster_label        : id and label at cluster level	
        supertype_id/supertype_label    : id and label at supertype level
        subclass_id/subclass_label	: id and label at subclass
        class_id/class_label	        : id and label at clas
        anatomical_annotation           : anatomical region that contribute to the most of cells in the cluster
        CCF_broad.freq/CCF_acronym.freq	: fraction of cells from regions in CCF_broad or CCF_acronym
        v3.size/v2.size/multiome.size   : number of cells or /nuclei from 10Xv3, 10Xv2, and multiome data
        cluster.markers                 : markers of the cluster
        merfish.markers                 : markers of the  cluster for merfish data
        max.region                      : The region that contribute to the most of the cells in given cluster
        max.region.ratio                : The fraction of cells coming from the max.region 
        F,M                             : sex distribution


