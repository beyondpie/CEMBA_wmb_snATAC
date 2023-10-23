cembav2env <- rlang::new_environment(
    data = list(
        # current L2 clustering updated at 220716
        L1CMeta = data.frame(
            L1Id = seq(3),
            L1Label = c("NonN", "GLUT", "GABA"),
            L1Reso = c(0.6, 1.0, 0.7)
        ),
        # * update L3 clustering with the selected resolution
        NonNL3CMeta = list(
            L2 = paste("NonN", seq(13), sep = "_"),
            L2Label = c("OGC", "ASC", "ASC",
                "OGC", "OPC", "MGL",
                "VLMC", "IOL", "VPIA",
                "PER", "TANY", "VEC",
                "ENTN"),
            resetL2 = c(
                S4Vectors::Pairs("COP", "IOL"),
                S4Vectors::Pairs("NFOL", "IOL"),
                S4Vectors::Pairs("MFOL", "OGC"),
                S4Vectors::Pairs("MOL", "OGC"),
                S4Vectors::Pairs("OPC", "OPC"),
                S4Vectors::Pairs("TANY", "EPEN"),
                S4Vectors::Pairs("EPMB", "EPEN"),
                S4Vectors::Pairs("HYPEN", "EPEN"),
                S4Vectors::Pairs("CHOR", "CHOR"),
                S4Vectors::Pairs("ASCG", "ASC"),
                S4Vectors::Pairs("ASCNT", "ASC"),
                S4Vectors::Pairs("ASCW", "ASC"),
                S4Vectors::Pairs("BERG", "BERG")
            ),
            L3Reso = c(1.0, 0.7, 1,
                0.3, 0.3, 0.3,
                0.5, 0.3, 0.3,
                0.5, 0.9, 0.5,
                1),
            L3nCluster = c(9, 8, 11,
                1, 3, 3,
                7, 3, 3,
                6, 11, 6,
                5),
            L3Label = list(
                NonN_1 = c(
                    S4Vectors::Pairs(3, "MFOL1"),
                    S4Vectors::Pairs(c(1:2, 4:9), paste0("MOL", 1:8))
                ),
                NonN_4 = S4Vectors::Pairs(1, "MOL9"),
                NonN_5 = S4Vectors::Pairs(1:3, paste0("OPC", 1:3)),
                NonN_6 = S4Vectors::Pairs(1:3, c("MGL1", "PVM", "MGL2")),
                NonN_7 = S4Vectors::Pairs(1:7, paste0("VLMC", 1:7)),
                NonN_8 = S4Vectors::Pairs(1:3, c("MFOL2", "COP", "NFOL")),
                NonN_9 = S4Vectors::Pairs(1:3, paste0("VPIA", 1:3)),
                NonN_10 = S4Vectors::Pairs(1:6, paste0("PER", 1:6)),
                NonN_12 = S4Vectors::Pairs(1:6, paste0("VEC", 1:6)),
                NonN_13 = S4Vectors::Pairs(1:5, paste0("ENTN", 1:5))
            )
        ),
        NonNL4CMetaExtra = data.frame(
            L4 = paste(paste("NonN", c(2, 3, 11), sep = "_", collapse = "-"),
                c(1:8), sep = "_"),
            L4Reso = c(0.3, 0.7, 0.9, 0.4, 1.0, 0.6, 0.5, 0.6),
            L4nCluster = c(1, 6, 10, 3, 11, 1, 2, 1)
        ),
        # * multiple-group combine result
        NonNL3CMetaExtra = list(
            L2 = c(paste("NonN", c(2, 3, 11), sep = "_", collapse = "-")),
            L2Label = c("ASC"),
            L3Reso = c(0.6),
            L3nCluster = c(8),
            resetL2 = S4Vectors::Pairs(
              "ASCG", "ASCG"
            ),
            L3Label = list(
                NonN_2.3.11 = c(
                    S4Vectors::Pairs(1, "ASCG"),
                    S4Vectors::Pairs(2, "ASCNT"),
                    S4Vectors::Pairs(3, "ASCW"),
                    S4Vectors::Pairs(4, "RGL"),
                    S4Vectors::Pairs(5, "EPEN"),
                    S4Vectors::Pairs(6, "ASCG"),
                    S4Vectors::Pairs(7, "BERG"),
                    S4Vectors::Pairs(8, "CHOR"))
            ),
            L4Label = list(
                NonN_2.3.11_1 = S4Vectors::Pairs(1, "ASCG1"),
                NonN_2.3.11_2 = S4Vectors::Pairs(1:6, paste0("ASCN", 1:6)),
                NonN_2.3.11_3 = S4Vectors::Pairs(1:10, paste0("ASCW", 1:10)),
                NonN_2.3.11_4 = S4Vectors::Pairs(1:3, c("RGSZ", "RGDG", "NIPC")),
                NonN_2.3.11_5 = c(
                    S4Vectors::Pairs(c(1,2,6,8,10), paste0("EPEN", 1:5)),
                    S4Vectors::Pairs(c(3,4,7,11), paste0("TANY", 1:4)),
                    S4Vectors::Pairs(5, "EPMB"),
                    S4Vectors::Pairs(9, "HYPEN")
                ),
                NonN_2.3.11_6 = S4Vectors::Pairs(1, "ASCG2"),
                NonN_2.3.11_7 = S4Vectors::Pairs(1:2, c("BERG1", "BERG2")),
                NonN_2.3.11_8 = S4Vectors::Pairs(1, "CHOR1")
            )
        ),

        GLUTL3CMeta = list(
            L2 = paste("GLUT",
                c(1, "2.7", 3, 4,
                    5, 6, "8.10", "9.12",
                    11, 13, 14, 15,
                    "16.25", 17, 18, 19,
                    20, 21, "22.24", 23,
                    26, 27, 28),
                sep = "_"),
            L2Label = c("GRANCB", "ITL23GL", "CTGL", "OLFGL",
                "GRC", "ITL4GL", "ITL6GL", "ITL5GL",
                "THGL", "PTGL", "CA1GL", "CA3GL",
                "ITHGL", "NPGL", "L6bGL", "PIRGL",
                "CLAGL", "GLUT_21", "AMY1GL", "GLUT_23",
                "GLUT_26", "GLUT_27", "AMY2GL"),
            resetL1 = NULL,
            resetL2 = c(
                S4Vectors::Pairs("OBGL", "OBGL"),
                S4Vectors::Pairs("OLFGL", "OFLGL"),
                S4Vectors::Pairs("PAGGL", "PAGGL"),
                S4Vectors::Pairs("THGL", "THGL"),
                S4Vectors::Pairs("GRANCB", "GRANCB"),
                S4Vectors::Pairs("GRANMY", "GRANMY")
            ),
            L3Reso = c(0.4, 1.0, 1.0, 1.0,
                0.8, 0.6, 0.5, 0.8,
                1.0, 0.7, 0.3, 0.4,
                0.4, 1.0, 1.0, 0.2,
                1.0, 1.0, 0.2, 0.2,
                0.3, 0.3, 0.0),
            L3nCluster = c(4, 9, 11, 17,
                7, 5, 7, 8,
                16, 11, 4, 8,
                5, 11, 10, 2,
                9, 6, 5, 2,
                1, 1, 1),
            L3Label = list(
              GLUT_1 = c(
              # GRANMB: Granule neurons, Cerebellum
              # GRANMY: Granule neurons, MY region
                S4Vectors::Pairs(1:3, paste0("GRANCB", 1:3)),
                S4Vectors::Pairs(4, "GRANMY1")
              ),
              GLUT_2.7 = S4Vectors::Pairs(1:12, paste0("ITL23GL", 1:12)),
              GLUT_3 = S4Vectors::Pairs(1:11, paste0("CTGL", 1:11)),
              GLUT_4 = c(
                # use OBGL 1,2,3 to label
                ## S4Vectors::Pairs(6, "OBGL1-2"),
                ## S4Vectors::Pairs(17, "OBGL3"),
                ## S4Vectors::Pairs(16, "OBGL4-5"),
                S4Vectors::Pairs(6, "OBGL1"),
                S4Vectors::Pairs(17, "OBGL2"),
                S4Vectors::Pairs(16, "OBGL3"),
                S4Vectors::Pairs(1:5, paste0("OLFGL", 1:5)),
                S4Vectors::Pairs(7:15, paste0("OLFGL", 6:14))
              ),
              GLUT_5 = S4Vectors::Pairs(1:7, paste0("GRC", 1:7)),
              GLUT_6 = S4Vectors::Pairs(1:5, paste0("ITL4GL", 1:5)),
              GLUT_8.10 = S4Vectors::Pairs(1:7, paste0("ITL6GL", 1:7)),
              GLUT_9.12 = S4Vectors::Pairs(1:8, paste0("ITL5GL", 1:8)),
              GLUT_11 = c(
                # PAG: Periaqueductal gray
                # APN: Anterior pretectal nucleus
                # MRN: Midbrain reticular nucleus
                ## GABA-class has 10 PAGGL 
                S4Vectors::Pairs(c(8,14), paste0("PAGGL", 11:12)),
                S4Vectors::Pairs(c(1:7, 9:13, 15:16), paste0("THGL", 1:14))
              ),
              GLUT_13 = S4Vectors::Pairs(1:11, paste0("PTGL", 1:11)),
              GLUT_14 = S4Vectors::Pairs(1:4, paste0("CA1GL", 1:4)),
              GLUT_15 = S4Vectors::Pairs(1:8, paste0("CA3GL", 1:8)),
              GLUT_16.25 = S4Vectors::Pairs(1:5, paste0("ITHGL", 1:5)),
              GLUT_17 = S4Vectors::Pairs(1:11, paste0("NPGL", 1:11)),
              GLUT_18 = S4Vectors::Pairs(1:10, paste0("L6bGL", 1:10)),
              GLUT_19 = S4Vectors::Pairs(1:2, paste0("PIRGL", 1:2)),
              GLUT_20 = S4Vectors::Pairs(1:9, paste0("CLAGL", 1:9)),
              GLUT_22.24 = S4Vectors::Pairs(1:5, paste0("AMY1GL", 1:5)),
              GLUT_28 = S4Vectors::Pairs(1, "AMY2GL1")
            )
        ),
        GABAL3CMeta = list(
            resetL1 = c(
              S4Vectors::Pairs("DGNBL", "GLUT"),
              S4Vectors::Pairs("MBRN", "CHOL;GLUT"),
              S4Vectors::Pairs("CHOLHA", "CHOL;GABA"),
              S4Vectors::Pairs("MBGL", "GLUT"),
              S4Vectors::Pairs("DOPMB", "DOP;GLUT"),
              S4Vectors::Pairs("SEROPAG", "SERO;GABA"),
              S4Vectors::Pairs("SEROPM", "SERO;GABA"),
              S4Vectors::Pairs("UnkIsoGA", "GABA"),
              S4Vectors::Pairs("ICGL", "GLUT"),
              S4Vectors::Pairs("PAGGL", "GLUT"),
              S4Vectors::Pairs("MYGL", "GLUT"),
              S4Vectors::Pairs("HYGL", "GLUT"),
              S4Vectors::Pairs("STRGL", "GLUT"),
              S4Vectors::Pairs("PALGL", "GLUT"),
              S4Vectors::Pairs("PEPT", "PEPT"),
              S4Vectors::Pairs("OREX", "OREX"),
              S4Vectors::Pairs("OXYT", "OXYT"),
              S4Vectors::Pairs("VASO", "VASO")
            ),
            L2 = paste("GABA", 1:33, sep = "_"),
            ## HBGA: hindbrain (pons and medulla) GABA
            L2Label = c("HBGA", "D12MSN", "HYGA", "PVGA", "MBGL",
                "SSTGA", "LSXGA", "PYRGA", "MBSC", "VIPGA",
                "STRGA", "LAMGA", "HYGL", "MSGA", "OBGA",
                "ICGL", "MBGA", "THGA", "OBNBL", "PONSGL",
                "PURKCB", "MXD", "STRGL", "DGNBL", "CHOLHA",
                "SEROPM", "PVGA", "STRGA", "DOPMB", "PVGA",
                "CRC", "GABA_32", "GABA_33"),
            resetL2 = c(
                S4Vectors::Pairs("D1MSN", "D1MSN"),
                S4Vectors::Pairs("D2MSN", "D2MSN"),
                S4Vectors::Pairs("OBDOP", "OBDOP"),
                S4Vectors::Pairs("OBNBL", "OBNBL"),
                S4Vectors::Pairs("STRGA", "STRGA"),
                S4Vectors::Pairs("MSGA", "MSGA"),
                S4Vectors::Pairs("CRC", "CRC"),
                S4Vectors::Pairs("MXD", "MXD"),
                S4Vectors::Pairs("PYRGA", "PYRGA"),
                S4Vectors::Pairs("MBGL", "MBGL"),
                S4Vectors::Pairs("MBRN", "MBRN"),
                S4Vectors::Pairs("MBSC", "MBSC"),
                S4Vectors::Pairs("THGA", "THGA"),
                S4Vectors::Pairs("PONSGA", "PONSGA"),
                S4Vectors::Pairs("MYGA", "MYGA"),
                S4Vectors::Pairs("MYGL", "MYGL"),
                S4Vectors::Pairs("HYGA", "HYGA"),
                S4Vectors::Pairs("MBGA", "MBGA"),
                S4Vectors::Pairs("MBSC", "MBSC"),
                S4Vectors::Pairs("MBRN", "MBRN"),
                S4Vectors::Pairs("MBGL", "MBGL"),
                S4Vectors::Pairs("PAGGL", "PAGGL"),
                S4Vectors::Pairs("HYGL", "HYGL"),
                S4Vectors::Pairs("STRGL", "STRGL"),
                S4Vectors::Pairs("PALGL", "PALGL"),
                S4Vectors::Pairs("SEROPAG", "SEROPAG"),
                S4Vectors::Pairs("SEROPM", "SEROPM"),
                S4Vectors::Pairs("PMCH", "PMCH"),
                S4Vectors::Pairs("PEPT", "PEPT"),
                S4Vectors::Pairs("OREX", "OREX"),
                S4Vectors::Pairs("OXYT", "OXYT"),
                S4Vectors::Pairs("VASO", "VASO"),
                S4Vectors::Pairs("CHOLHA", "CHOLHA"),
                S4Vectors::Pairs("RMSTTH", "RMSTTH"),
                S4Vectors::Pairs("PURKCB", "PURKCB"),
                S4Vectors::Pairs("PURKMY", "PURKMY"),
                S4Vectors::Pairs("KF", "KF"),
                S4Vectors::Pairs("UnkIsoGA", "UnkIsoGA")
            ),
            L3Reso = c(1.0, 1.0, 1.0, 0.9, 1.0,
                0.6, 1.0, 1.0, 1.0, 0.3,
                1.0, 0.6, 1.0, 1.0, 0.3,
                1.0, 0.6, 1.0, 0.4, 0.2,
                0.8, 0.5, 0.3, 0.5, 0.5,
                0.9, 0.2, 0.5, 0.6, 0.3,
                0.6, 0.0, 0.0),
            L3nCluster = c(30, 12, 29, 11, 30,
                10, 21, 24, 28, 5,
                11, 6, 25, 20, 4,
                16, 12, 9, 4, 1,
                8, 2, 5, 2, 4,
                15, 1, 1, 7, 2,
                1, 1, 1),
            L3Label = list(
              GABA_1 = c(
                S4Vectors::Pairs(c(1,6, 27), paste0("MYGA", 1:3)),
                S4Vectors::Pairs(c(2,4,5,7,9,10, 17,19,24), paste0("MBGL", 28:36)),
                # SEROPAG: Serotonergic neurons, hindbrain PAG
                S4Vectors::Pairs(18, "SEROPAG1"),
                S4Vectors::Pairs(c(16,20), paste0("RMSTTH", 1:2)),
               # SEROPM: Serotonergic neurons in Pons, Medulla
                S4Vectors::Pairs(c(3,8,23,25), paste0("SEROPM", 6:9))
              ),
              # NOTE: CNUGA locates in GABA_2_1, GABA_2_7 mainly
              # CNUGA closed to D2MSN.
              GABA_2 = c(
                S4Vectors::Pairs(c(1:3, 6:7), paste0("D2MSN", 1:5)),
                S4Vectors::Pairs(c(4,5, 8:12), paste0("D1MSN", 1:7))
                ## FROM UMAP, cluster 9 did not seperate far away
                ## S4Vectors::Pairs(9, "STRGA13")
              ),
              GABA_3 = c(
                S4Vectors::Pairs(c(2,21), c("PEPT1", "PEPT2")),
                S4Vectors::Pairs(28, "OREX1"),
                S4Vectors::Pairs(c(24,27), paste0("OXYT", 1:2)),
                S4Vectors::Pairs(1, "VASO1"),
                S4Vectors::Pairs(c(3:17,20, 23, 26, 29), paste0("HYGA", 1:19)),
                S4Vectors::Pairs(c(19,25), paste0("RMSTTH", 3:4)),
                S4Vectors::Pairs(18, "MSGA21"),
                S4Vectors::Pairs(22, "MSGA22")
              ),
              GABA_4 = S4Vectors::Pairs(1:11, paste0("PVGA", 1:11)),
              GABA_27 = S4Vectors::Pairs(1, "PVGA12"),
              GABA_30 = S4Vectors::Pairs(1:2, paste0("PVGA", 13:14)),
              ## PAG: periaqueductal gray,
              ## MRN: Midbrain reticular nucleus
              GABA_5 = c(
                S4Vectors::Pairs(c(1:14), paste0("MBGL", 1:14)),
                S4Vectors::Pairs(c(16:22), paste0("MBGL", 15:21)),
                S4Vectors::Pairs(c(25:30), paste0("MBGL", 22:27)),
                ## MBRN: cholinergic neurons Midbrain red nucleus
                S4Vectors::Pairs(c(15,23,24), paste0("MBRN", 1:3))
              ),
              GABA_9 = c(
                S4Vectors::Pairs(c(1,3,4,5,7,8,10,11,12,14,15,16,20,22,23),
                  paste0("MBGA", 1:15)),
                # MBSC: superior colliculus, midbrain
                S4Vectors::Pairs(c(2,6, 9,17,18,19,24,25,28),
                  paste0("MBSC", 1:9)),
                S4Vectors::Pairs(c(13,21), paste0("RMSTTH", 5:6)),
                S4Vectors::Pairs(26, "PONSGA1"),
                S4Vectors::Pairs(27, "PMCH1")
              ),
              GABA_29 = c(
                # DOPMB: Dopaminergic neurons, ventral midbrain (SNc, VTA)
                S4Vectors::Pairs(c(1:3, 6,7), paste0("DOPMB", 1:5)),
                # Koelliker-Fuse subnucleus
                S4Vectors::Pairs(4:5, paste0("KF", 1:2))
              ),
              GABA_29 = S4Vectors::Pairs(1:7, paste0("DOPMB", 1:7)),
              GABA_17 = c(
                S4Vectors::Pairs(c(1:10, 12), paste0("MBGA", 16:26)),
                # UnkIsoCGA: unknown Isocortext GA
                S4Vectors::Pairs(11, "UnkIsoGA1")
              ),
              # IC: Inferior colliculus, midbrain
              GABA_16 = S4Vectors::Pairs(1:16, paste0("ICGL", 1:16)),
              GABA_13 = c(
                S4Vectors::Pairs(c(1,3, 5,6,9,10,11,14,18,21),
                  paste0("PAGGL", 1:10)),
                S4Vectors::Pairs(c(2,7,15,16,19,20),
                  paste0("MYGL", 1:6)),
                S4Vectors::Pairs(c(4,17,23,25), paste0("HYGL", 1:4)),
                S4Vectors::Pairs(8, "STRGL1"),
                ## PALGL" PAL-region majored GLUT
                S4Vectors::Pairs(c(12, 22, 24), paste0("PALGL", 1:3))
              ),
              GABA_6 = S4Vectors::Pairs(1:10, paste0("SSTGA", 1:10)),
              GABA_28 = S4Vectors::Pairs(1, "STRGA12"),
              GABA_7 = S4Vectors::Pairs(1:21, paste0("LSXGA", 1:21)),
              GABA_10 = S4Vectors::Pairs(1:5, paste0("VIPGA", 1:5)),
              GABA_12 = S4Vectors::Pairs(1:6, paste0("LAMGA", 1:6)),
              # and have a markder Grp
              GABA_14 = c(
                S4Vectors::Pairs(1:5, paste0("MSGA", 1:5)),
                # mapped to PYRGA9
                S4Vectors::Pairs(6, "PYRGA24"),
                S4Vectors::Pairs(7:20, paste0("MSGA", 6:19))
              ),
              GABA_15 = S4Vectors::Pairs(1:4, paste0("OBGA", 1:4)),
              GABA_19 = c(
                S4Vectors::Pairs(c(1,3), paste0("OBNBL", 1:2)),
                S4Vectors::Pairs(c(2,4), paste0("OBDOP", 1:2))
              ),
              GABA_11 = S4Vectors::Pairs(1:11, paste0("STRGA", 1:11)),
              GABA_22 = S4Vectors::Pairs(1, "MXD1"),
              GABA_24 = S4Vectors::Pairs(1:2, paste0("DGNBL", 1:2)),
              GABA_31 = S4Vectors::Pairs(1, "CRC1"),
              # this aligned with AMY region PYRGA family
              GABA_8 = c(
                S4Vectors::Pairs(5, "D2MSN6"),
                S4Vectors::Pairs(1:4, paste0("PYRGA", 1:4)),
                S4Vectors::Pairs(6:24, paste0("PYRGA", 5:23))
              ),
              GABA_18 = S4Vectors::Pairs(1:9, paste0("THGA", 1:9)),
              # CHOLHA: cholinergic neurons, habenula
              GABA_25 = S4Vectors::Pairs(1:4, paste0("CHOLHA", 1:4)),
              GABA_21 = c(
                S4Vectors::Pairs(c(1:5, 7:8), paste0("PURKCB", 1:7)),
                S4Vectors::Pairs(6, "PURKMY1")
              ),
              GABA_20 = S4Vectors::Pairs(1, "PONSGL1"),
              # SEROPM: Serotonergic neurons in Pons, Medulla
              # left GABA_26 cluster: 6-15, currently label their L2
              # as SEROPM
              GABA_26 = S4Vectors::Pairs(1:5, paste0("SEROPM", 1:5))
            )
        ),
        classMap = c(
          "NN" = "NN", "Glut"= "Glut", "Gaba" = "GABA",
          "Dopa" = "Dopa", "Nora" = "Nora", "Sero" = "Sero",
          "Gaba;Glut"  = "GABA-Glut", "Glut;Gaba" = "GABA-Glut",
          "Gaba;Chol" = "GABA-Chol",
          "Gaba;Dopa" = "GABA-Dopa", "Gaba;Gly" = "GABA-Gly",
          "Glut;Chol" = "Glut-Chol", "Glut;Dopa" = "Glut-Dopa",
          "Gaba;Hist" = "GABA-Hist"),
        anatomicalMap = list(
          Telencephalon = c("Isocortex", "HPF", "OLF", "AMY", "STR", "PAL"),
          Diencephalon = c("TH", "HY"),
          Midbrain = "MB",
          Hindbrain = c("Pons", "MY"),
          Cerebellum = "CB"),
        # * cembav1 meta
        cembav1MetaFile = file.path(
            here::here(),
            "meta", "rs1cemba.full.cellMeta.taxonomyLable.tsv"),
        # raw snap with no quality control and doublets.
        qcTaijiDir = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "qc.Taiji"),
        snapListBeforeQC = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "snap.mapping.meta.before.QC.rds"),
        snMetaBeforeQC = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "snMetaBeforeQC.v1.rds"),
        doubletStatDir = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "doubletStat"),
        doubletDir = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "doublets"),
        # current L2 clustering result dir
        snapL2Dir = file.path(here::here(),
          "02.annotation", "workflow/RANNL2/out"),
        ## current L3 clustering result dir
        snapL3Dir = file.path(here::here(), "02.annotation", "fossa/L3_Result"),
        snapAllFile = file.path(here::here(), "02.annotation", "fossa",
            "snap", "snapWithL3Extra.rds"),
        brainRegionFile = file.path(here::here(), "meta",
          "BrainRegion.Metadata.txt"),
        dissect2TimeFile = file.path(here::here(), "meta", "dissect2time.csv"),
        # * annotation-related files
        L3ClassInfoFile = file.path(here::here(),
          "meta", "L3.clas.info.csv"),
        # * current doublet files
        pmatDoubletFile = file.path(here::here(), "supple.00.doublets",
          "CEMBAv1.pmat.refineDoublets.txt"),
        gmatDoubletFile = file.path(here::here(), "supple.00.doublets",
            "CEMBA.gmat.refineDoublets.txt"),
        gmatRawDoubletFile = file.path(here::here(), "supple.00.doublets",
            "CEMBA.gmat.raw.fitDoublets.txt"),
        # * current meta all in one dir
        metaDir = file.path(here::here(), "supple.02.annotation.all.in.one.meta",
          "data"),

        ## - meta with Integration v1, where GLUT/GABA
        ## had inconsistent labels at that time
        ## on main class level
        ## allmeta.v6.File = file.path(here::here(),
        ##   "supple.02.annotation.all.in.one.meta",
        ##   "data", "cell.mergeASC", "mba.whole.cell.meta.v6.4.rds"),
        
        ## - add allen annotation based on Integration on L3-level
        ## (mainclass-based) for RNAseq and Integration of snMethy
        ## - update biorep, which has bug on v7.1 for 15B
        ## allmetaFile = file.path(here::here(),
        ##   "supple.02.annotation.all.in.one.meta",
        ##   "mba.whole.cell.meta.v7.2.rds"),

        ## - add L2 annot as L2Annot2
        ## allmetaFile = file.path(here::here(),
        ##   "supple.02.annotation.all.in.one.meta",
        ##   "mba.whole.cell.meta.v7.3.rds"),

        ## - update GABA_1_29 class as Gaba;Hist
        ## - add columns: classv2 and classv2Color
        ## allmetaFile = file.path(here::here(),
        ##   "supple.02.annotation.all.in.one.meta",
        ##   "mba.whole.cell.meta.v7.4.rds"),
        ## - Change PVGA2 -> PVGL, PVGA1 -> PVGA
        ## allmetaFile = file.path(here::here(),
        ##   "supple.02.annotation.all.in.one.meta",
        ##   "mba.whole.cell.meta.v7.5.rds"),
        ## add doublet score
        allmetaFile = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "mba.whole.cell.meta.v8.1.rds"),
        
        L3MetaOldFile = file.path(here::here(), "meta", "L3Meta.csv"),
        L3MetaFile = file.path(here::here(),
          "meta", "L3Meta.v7.1.txt"),
        L3OrderFile = file.path(here::here(), "meta",
          "L3toAllen.hc.order.txt"),
        # * umap file
        umapFile = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "cell.umap", "cemba.all.umap.matrix.rds"),
        # * current gmat
        gmatDir = file.path(here::here(), "supple.02.annotation.all.in.one.meta",
          "data", "cell.gmat"),
        gmatFileSuffix = ".dp-3e+06.gmat.rds",
        gmatvM23Dir = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "cell.gmat.vM23"),
        mm.gencode.vM23.gtf = file.path(here::here(),
          "meta", "modified_gencode.vM23.primary_assembly.annotation.gtf"),
        mm.gene.tss.updn1k.vM23.bed = file.path(here::here(),
          "meta", "gencode.vM23.gene.tssUpDn1k.bed"),
        ens2symbolFile = file.path(here::here(), "meta/ensemble2symbol_vM23.tsv"),
        collectedGeneMarkerFile = file.path(here::here(), "meta/geneMarkerV2.csv"),
        # * integration
        mc2atacFile = file.path(here::here(), "meta/mc2atac.Int.csv"),
        # * original region-specific annoation dir
        regionClusteringDir = file.path(here::here(), "01.cluster",
          "regionAnnotationFromYang"),
        L2Annot2File = file.path(here::here(),
          "meta", "L2Annot2.txt"),
        L2HCFile = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "L2Annot2.hc.pearson.logcpm.cov0.5.rds"),
        cellClusterMetav1 = file.path(here::here(),
          "meta", "cellcluster.annot.CEMBAv1.txt"),
        # * peak-related information
        peakDir = file.path(here::here(), "supple.07.peakcalling.allinone"),
        peakBedFile = file.path(here::here(), "supple.07.peakcalling.allinone",
          "mba.whole.fitPeakModel.q0.01.union.bed"
        ),
        cCREBedFile = file.path(here::here(), "supple.07.peakcalling.allinone",
          "whole.mouse.brain.cCREs.bed"
        ),
        ovlpDHSPeakBedFile = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.q0.01.OvlpDHS.bed"
        ),
        nonOvlpDHSPeakBedFile = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.q0.01.nonOvlpDHS.bed"
        ),
        shufflePeakBedFile = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.shuffle.removeovlp.bed"
        ),
        is.ovlpDHS.File = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "is.ovlpDHS.allPeaks.csv"
          ),
        ## subclassPmatCPMIntv1File = file.path(here::here(),
        ##   "supple.07.peakcalling.allinone",
        ##   "mba.whole.subclass.pmat.sc.filter.cpm.Intv1.rds"),
        subclassPmatCPMIntv2File = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.subclass.pmat.sc.filter.cpm.IntL3v2.rds"),
        subclassPmatCPMIntv2File2 = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.subclass.pmat.sc.filter.cpm.IntL3v2.peaknm.rds"),
        subclassPmatCountIntv2File = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.subclass.pmat.sc.filter.totalcount.rds"
        ),
        subclassPmatBinaryIntv2File = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.subclass.pmat.binary.rds"
        ),
        allenClassAndRegionIntv2 = file.path(here::here(),
          "meta", "allen2region.60pcnt.csv"),
        allen2class2mainclassFile = file.path(here::here(),
          "meta", "allenAnnot2class2mainclass.csv"
          ),
        homerPeakAnnotSumFile = file.path(
          here::here(), "supple.07.peakcalling.allinone",
          "homer.annot.scfilter.peak.txt"
        ),
        IntL3MetaFile = file.path(here::here(), "meta", "IntL3SumMeta.csv"),
        celltypePmatCPMFile = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.celltype.pmat.sc.filter.cpm.rds"),
        celltypePmatBinaryFile = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.celltype.pmat.binary.rds"),
        celltypePmatCountFile = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.celltype.pmat.sc.filter.cnt.rds"),
        L2PmatCountFile = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mba.whole.L2.pmat.sc.filter.cnt.rds"),
        L3PeakDir = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "L3PeakFinal3"),
        embedDir = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "cell.embed"),
        embedL1File = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "cell.embed",
          "cemba.all.embed.matrix.L1.rds"),
        embedL2File = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "cell.embed",
          "cemba.all.embed.matrix.L2.rds"),
        embedL3File = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "data", "cell.embed",
          "cemba.all.embed.matrix.L3.rds"),
        annotFile = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "atac.int.summary.v1.tsv"),
        allenGeneFile = file.path(here::here(), "meta",
          "ensemble.genesymbol.allengenesymbol.csv"),
        bigwigL3Dir = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "celltype_bw"),
        allenCell2Region2L2File = file.path(here::here(), "03.integrate",
          "allen.data", "cell2region2l2.allen.10xv3.male.all.rds"),
        allGeneSymbolSharedAtacAllenFile = file.path(here::here(), "meta",
          "gene.symbol.shared.atac.allen.csv"),
        allenClusteringMetaFile = file.path(here::here(), "03.integrate",
          "allen.data", "allen.clustering.meta.szu.csv"),
        allenRegionToCembaMajorRegionFile = file.path(here::here(), "meta",
          "allen.region.to.main.region.txt"),
        allenSubclassMetaFile = file.path(here::here(), "meta",
          "allen.subclass.meta.csv"),
        allen.l2.10xv3.nn.ds200.seurat = file.path(here::here(), "03.integrate",
          "before_integration", "allen.subclass.nn.ds.200.seurat.rds"),
        allen.l2.10xv3.neuron.ds200.seurat = file.path(
          here::here(), "03.integrate",
          "before_integration", "allen.subclass.neuron.ds.200.seurat.rds"),
        allenAnnoctConcatList = file.path(here::here(),
          "supple.02.annotation.all.in.one.meta",
          "mba.whole.allenAnnotConcat.list"),
        allen.l2.10xv3.ds200.cpm = file.path(here::here(),
          "04.cCREgene", "allen.l2.cpm.ds200.rds"),
        # 215 subclasses in total
        subclass.order.hc.Intv2 = file.path(here::here(),
          "meta", "subclass.order.hc.Intv2.csv"),
        # mCG signals
        mCG.L3.mean.File = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mCGmat.mean.rds"),
        mCG.all.File = file.path(here::here(),
          "supple.07.peakcalling.allinone",
          "mCG.allpeak.rds"),
        # cCRE-gene pairs
        all.pdc.info = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "mba.whole.AllenAnnotConcat.merge.pdc.all"),
        all.pdc.pairs.only = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "mba.whole.AllenAnnotConcat.merge.pdc.pair.all"),
        all.pdc.dist = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "mba.whole.AllenAnnotConcat.merge.pdc.dist.all"),
        pdc.cor.spearman.real.file = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "AllenAnnotConcat.cor.spearman.real.csv"),
        pdc.cor.spearman.rdmshuf.file = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "AllenAnnotConcat.cor.spearman.rdm.shuf.csv"),
        pdc.cor.spearman.fdr.tsv = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "AllenAnnotConcat.pdc.cor.spearman.fdr.alignv1.tsv"),
        pdc.cor.pearson.real.file = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "AllenAnnotConcat.cor.pearson.real.csv"),
        pdc.cor.pearson.rdmshuf.file = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "AllenAnnotConcat.cor.pearson.rdm.shuf.csv"),
        pdc.cor.pearson.fdr.tsv = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "AllenAnnotConcat.pdc.cor.pearson.fdr.alignv1.tsv"),
        ppdc.pearson.tsv = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "mba.whole.AllenAnnotConcat.pearson.pos.pdc.alignv1.tsv"),
        npdc.pearson.tsv = file.path(here::here(),
          "04.cCREgene", "out", "AllenAnnotConcat",
          "mba.whole.AllenAnnotConcat.pearson.neg.pdc.alignv1.tsv"),
        GRNmatsFile = file.path(here::here(),
          "05.GRN", "GRNmats.topk5000.nlogp2.v1.rds")
    ),
    parent = rlang::empty_env()
)
