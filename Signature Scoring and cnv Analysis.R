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
#####V## GO KEGG###
suppressMessages({
  library(Seurat)
  #library(corrplot)
  library(tidyverse)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  #library(org.Mm.eg.db)
})

Idents(scRNA1)<-"celltype"
scRNA2<-subset(scRNA1, idents=c('FPR_Mo'))
cluster1.markers <- FindMarkers(scRNA2, ident.1 = 'Tumor',ident.2 = 'HI', verbose = TRUE, logfc.threshold=0.05)
sig_dge.celltype <- subset(cluster1.markers, p_val_adj<0.05&avg_log2FC>0.5)
ego_ALL <- enrichGO(gene          = row.names(sig_dge.celltype),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)   
ego_CC <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(sig_dge.celltype),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP,showCategory = c('neutrophil chemotaxis','positive regulation of inflammatory response','antigen processing and presentation of exogenous antigen','neutrophil migration','leukocyte migration','myeloid leukocyte migration','granulocyte migration')) + ggtitle("barplot for Biological process")
library(ggplot2)
selected_terms <- c('viral process', 'viral life cycle', 
                    'antigen processing and presentation of exogenous antigen', 
                    'regulation of viral process', 'antigen processing and presentation'
                    )
plot_data <- ego_BP@result %>% 
  filter(Description %in% selected_terms) %>% 
  arrange(desc(Count)) %>% 
  mutate(Description = factor(Description, levels = Description))
ggplot(plot_data, aes(x = Count, y = reorder(Description, Count))) +
  geom_segment(aes(x = 0, xend = Count, yend = Description),
               color = "#B383B9", linewidth = 1.2) +
  geom_point(size = 10, color = "#B383B9") +
  geom_text(aes(label = Count), 
            nudge_x = 0, 
            color = "white",  
            size = 4.5, 
            fontface = "bold") +  
  labs(title = "Biological Process Lollipop Plot",
       x = "Gene Count", y = "") +
  theme_minimal() +
  theme(
    axis.line  = element_line(color = "black", linewidth = 0.6), 
    axis.ticks  = element_line(color = "black"), 
    panel.grid  = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.x  = element_text(size = 11, face = "plain"),
    axis.title.x  = element_text(size = 13, margin = margin(t=10)),
    plot.title  = element_text(size = 16, hjust = 0.5, face = "bold.italic") 
  )

#######PROGENy#######
library(progeny)
library(tidyr)
library(tibble)
Idents(scRNA2) <- 'celltype'
CellsClusters <- data.frame(Cell = names(Idents(scRNA2)), 
                            cellType = as.character(Idents(scRNA2)),
                            stringsAsFactors = FALSE)
scRNA2<- progeny(scRNA2, scale=FALSE, organism="Human", top=500, perm=1, assay_name = "SCT",
                 return_assay = TRUE)
scRNA2@assays$progeny
scRNA2<- Seurat::ScaleData(scRNA2, assay = "progeny") 
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(scRNA2, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
dim(progeny_scores_df)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, cellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
dim(summarized_progeny_scores)
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
paletteLength = 100
myColor = colorRampPalette(c("#4991C1", "white","#FF0000"))(paletteLength)
zzm8colors<-c('#C6307C',
              '#4991C1',
              '#D0AFC4',
              '#89558D',
              '#AFC2D9',
              '#435B95',
              '#79B99D','#E9E55A','#88558D','#C6367A',
              '#B3B5B1','#8BAACD',
              '#8E9296')

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

progeny_hmap = pheatmap(t(summarized_progeny_scores_df),fontsize=12, scale="column",cluster_rows=T,cluster_cols=F,
                        fontsize_row = 10, 
                        color=myColor)

















