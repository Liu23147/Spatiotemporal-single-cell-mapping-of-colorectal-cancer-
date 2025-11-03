library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
gene_list <- list(
  MyGenes = c("NR1H3","SIRPB1","ABCA7","RACK1","SPON2","VAV3","MERTK","CEBPE","COLEC10","CEACAM4","LMAN2","RAB31","CORO1A","TREX1","LYST","CD300A","TUSC2","XKR4","ARAP1","CLCN2","CLCN3","CLN3","PLD4","SPACA3","CNN2","ADORA1","ADORA2A","CRP","SIRPA","CRYBA1","CSK","CD300LF","TAFA4","CYBA","AGER","DNM2","DOCK1","DOCK2","ABCA1","PLPP4","ANO6","AHSG","AIF1","ELANE","PIKFYVE","UNC13D","F2RL1","FCER1G","FCER2","FCGR1A","FCGR1BP","FCGR2A","FCGR2B","FCN1","FCN2","FGR","RAB11FIP2","CD93","MESD","SYT11","JMJD6","PIP5K1C","FPR2","CORO1C","SH3BP1","ALOX15","ABL1","FYN","LETMD1","C1orf43","APPL1","GAS6","STAP1","GATA2","VPS33B","HAVCR1","EIF2AK1","XKR6","PYCARD","GSN","ANXA1","HCK","ANXA3","NCKAP1L","ANXA11","HMGB1","APOA1","APOA2","RAB7B","ICAM3","IRF8","XKR7","IFNG","TICAM2","IL2RB","IL2RG","IL15","IL15RA","ITGA2","ITGAL","ITGAM","ITGAV","ITGB1","ITGB2","ITGB3","PEAR1","LBP","LDLR","LEP","LEPR","LIMK1","MYO18A","LRP1","LYN","MIR17","MIR181B1","MIR183","MIR20A","MBL2","MFGE8","MSR1","MST1R","MYD88","MYH9","MYO7A","NCF2","NCF4","NOS2","P2RX7","P2RY6","PAK1","BIN2","GULP1","RAB14","PECAM1","CFP","PIK3CA","PLA2G5","PLCG2","PLD2","PLSCR1","TREM2","RAB39A","XKR8","APPL2","ARL8B","SIRPG","PIP4P2","LYAR","RAB20","SLC48A1","PRKACA","AXL","PRKCD","PRKCE","PRKCG","PRTN3","AZU1","CDC42SE1","CDC42SE2","CAMK1D","PTK2","ADGRB1","PTPRC","PTPRJ","PTX3","RAB5A","RAB27A","RAC1","RAP1A","RAP1GAP","RARA","BCR","CCL2","VIPAS39","ELMO2","MYO1G","SFTPD","ATG3","CLEC7A","SLAMF1","SFTPA1","SLC11A1","SOD1","SRC","STAT3","VAMP7","SYK","BTK","TGFB1","TGM2","THBS1","ICAM5","TLR2","TLR4","C2","C3","C4A","C4B","C4BPA","C4BPB","TUB","TULP1","CCR2","TYRO3","TYROBP","VAV1","VAV2","YES1","RAB7A","COLEC11","ELMO3","SPG11","COLEC12","CALR","DYSF","RAB34","PIP5K1A","PLA2G6","SRPX","BLTP1","TMEM175","MEGF10","FCN3","BECN1","MARCO","SNX3","SPHK1","SYT7","FCGR2C","TIMD4","CD14","ADIPOQ","ARHGAP12","ATG5","CD36","SCARB1","CD47","RUBCN","TM9SF4","ELMO1","CD302","ARHGAP25","CDC42")
)
scRNA3 <- AddModuleScore(object = scRNA3, features = gene_list, ctrl = 100, name = "Phagocytic_score")
meta_data <- scRNA3@meta.data
zzm8colors<-c('#C6307C',
              '#4991C1',
              '#D0AFC4',
              '#89558D',
              '#AFC2D9',
              '#435B95',
              '#79B99D')
meta_data$celltype <- factor(meta_data$celltype)
lot(meta_data, aes(x = celltype, y = Signature Scoring1)) +  
  geom_boxplot(aes(fill = celltype), outlier.shape = NA) +  
  labs(title = "Signature Scoring",
       x = "celltype",
       y = "Signature Scoring") +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"  
  )+ 
  scale_fill_manual(values =zzm8colors)
library(infercnv)
dfcount = as.data.frame(scRNA1@assays$SCT@counts)
groupinfo= data.frame(cellId = colnames(dfcount))
cellInfo <- pbmc_subset@meta.data[,c("tissue","celltype")]
colnames(cellInfo) <- c('tissue', 'celltype')
head(cellInfo)
head(cellInfo)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dfcount,
                                    annotations_file=cellInfo,
                                    delim="\t",
                                    gene_order_file="https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt",
                                    ref_group_names='HI')
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="output_dir", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)
