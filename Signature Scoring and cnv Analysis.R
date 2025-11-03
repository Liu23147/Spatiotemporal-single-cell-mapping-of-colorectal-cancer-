library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
gene_list <- list(
  HYPOXIA_gene = c(NDRG1,LDHA,PGK1,TPI1,PGAM1,PFKP,ENO1,SLC2A1,PNP,BNIP3,PSRC1,SLC16A1,DDIT4,ADM,HK2,SEC61G,AK4,CA9,P4HA1,ACOT7,GPI,TUBB6,TUBA1C,CDKN3,CTSV,VEGFA,TUBA1A,LRRC42,PSMA7,GAPDH,CHCHD2,YKT6,MIF,MAP7D1,MRPL15,MRPL13,MCTS1,UTP11,KIF4A,HILPDA,MRGBP,KIF20A,MRPS17,ESRP1,SHCBP1,SLC25A32,CORO1C,ANLN,MAD2L2,ANKRD37)
,M2_gene =(MMP9, MS4A6A, CCL18, AIF1, NCF2, CD4, ACP5, CLEC4A, CCL13, LY86, SLC15A3, CLEC10A, CCL23, HLA-DQA1, RNASE6, ADAMDEC1, HCK, NPL, TREM2, IRF8, SIGLEC1, CCL4, SAMSN1, CLEC7A, TLR2, P2RY13, CCL8, CD86, CD180, CD68, CD209, C3AR1, GPR183, CCL14, DPEP2, FPR3, CD37, LST1, CLIC2, HRH1, EGR2, CHI3L1, C1orf54, CFP, C5AR1, MNDA, PTGER2, CD300A, PIK3IP1, IGSF6, RENBP, PLA2G7, FCGR2B, SLAMF8, RASGRP3, QPCT, LILRB2, ATP8B4, PTGIR, FRMD4A, TLR8, KYNU, CHST15),M1_gene =(CXCL9, CCL19, CXCL10, IDO1, EBI3, TNFAIP6, CCR7, LAMP3, CCL5, CD40, IRF8, CD38, CCL4, APOL3, CXCL11, ADAMDEC1, HCK, BCL2A1, RSAD2, CYP27B1, KYNU, SIGLEC1, SLAMF1, SLC2A6, IFIH4, SAMSN1, IL2RB, TLR2, CCL8, APOBEC3A, AQP9, C3AR1, HLA-DQA1, MNDA, LILRA2, SLC15A3, CD80, ADGRE1, CHI3L1, RASSF4, IL4R, AIM2, PLA1A, MMP9, TLR8, HESX1, SOCS1, LST1, NOD2, IL7R, FLVCR2, BIRC3, PTGIR, CHI3L2, PLA2G7, GPR183, AIF1, CD4, HLA-DOB, MSC, ADGRE2, CCR5, CLEC2D, TNFSF4, CD86, SPIB, LAG3, DHX58, CLIC2, ACHE, CCL18, LY86, CCL14, TNIP3, C1orf54, MMP25, ACP5, APOBEC3G),CYTOLYTIC_gene = c(GNLY,RNF19B,SERPINB1,GZMH,GZMB,KIF5B,LAMP1,NKG7,ARL8B,PRF1,SRGN,STXBP2,TIAL1,CALR),ANGIOGENESIS_gene =(ACVRL1,AGGF1,AMOT,ANG,ANGPTL3,ANGPTL4,ATP5IF1,BTG1,C1GALT1,CANX,CDH13,CHRNA7,COL4A2,COL4A3,EGF,EMCN,EPGN,ERAP1,FOXO4,HTATIP2,IL17F,IL18,CXCL8,MYH9,NCL,NF1,NOTCH4,NPPB,NPR1,PF4,PLG,PML,PROK2,RHOB,RNH1,ROBO4,RUNX1,SCG2,SERPINF1,SHH,SPHK1,SPINK5,STAB1,TGFB2,THY1,TNFSF12,TNNI3,VEGFA))
scRNA1 <- AddModuleScore(object = scRNA1, features = gene_list, ctrl = 100, name = "score")
meta_data <- scRNA1@meta.data
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
