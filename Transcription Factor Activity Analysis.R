
library(nichenetr)
library(Seurat)
library(tidyverse)
library(circlize)
library(nichenetr)
library(viridis)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(SCENIC)
scRNA_subset=readRDS('merged_seurat.rds')

exprMat <- as.matrix(scRNA_subset@assays$SCT@data)#表达矩阵

cellInfo <- scRNA_subset@meta.data[,c("celltype","nCount_RNA","nFeature_RNA")]
colnames(cellInfo) <- c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", nCores=20)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
p<-pheatmap(regulonActivity_byCellType_Scaled,treeheight_row=15, treeheight_col=5, border_color=NA,color=colorRampPalette(c("#E7FBB9","white","#A73E28"))(100), scale = "row")
