####nichenetr####
library(pheatmap)
library(nichenetr)
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj =scRNA10, 
  receiver = c('CXCL13-Th',  'FCGR3A-NK',  'GZMK-Teff',   'HSPA1A-T',     'IFNG-T',    'IL17-Th',    'Naive-T','Treg',    'XCL1-NK'), 
    condition_colname = "tissue", condition_oi = c('ANT','Tumor'), condition_reference = 'HI',geneset='up',assay_oi='SCT',
  sender = c('FPR_Mo'),
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks,top_n_ligands=20)
pp=nichenet_output[["ligand_target_heatmap"]][["data"]]
result <- pp[pp$y == "ANXA1", ]
result <-subset(result,score>0.008)
cols_to_keep <- sapply(result, function(col) !all(col == 0))
df_cleaned <- result[, cols_to_keep]

p = nichenet_output$ligand_activity_target_heatmap
pheatmap(result)
library(pheatmap)


matrix_data <- result %>% 
  pivot_wider(names_from = x, values_from = score) %>% 
  column_to_rownames("y") %>% 
  as.matrix() 


pheatmap(matrix_data,
         color = colorRampPalette(c("#E1BEE7", "#4A148C"))(50),
         cluster_rows = FALSE,  
         cluster_cols = T,  
         display_numbers = F, 
         number_format = "%.2e", 
         fontsize_row = 12,     
         fontsize_col = 10,    
         angle_col = 45,        
         main = "Top genes of FPR_Mo cells")
nichenet_output %>% names()
nichenet_output$ligand_activities
nichenet_output$top_ligands
p = nichenet_output$ligand_target_heatmap
p = nichenet_output$ligand_activity_target_heatmap
p = nichenet_output$ligand_receptor_heatmap
p = nichenet_output$ligand_receptor_df
#######CellChat######################
library(CellChat)
Idents(scRNA1)<- 'tissue'
ANT<-subset(scRNA1, idents=c("ANT"))
Tumor<-subset(scRNA1, idents=c("Tumor"))
Idents(ANT)<- 'celltype'
Idents(Tumor)<- 'celltype'
data_obj <- ANT
data_obj$celltype <- Idents(data_obj)
data_obj
cellchat <- createCellChat(object = data_obj, meta = data_obj@meta.data, group.by = "celltype") 
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 
groupSize
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchatANT <- cellchat
data_obj <- Tumor
data_obj$celltype <- Idents(data_obj)
data_obj
cellchat <- createCellChat(object = data_obj, meta = data_obj@meta.data, group.by = "celltype") 
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 
groupSize
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchatTumor <- cellchat

object.list <- list(STK = cellchatANT,DMSO = cellchatTumor)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

gg1 <- netVisual_heatmap(cellchat)
## Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
## Do heatmap based on a merged object
gg1 + gg2

netVisual_bubble(cellchat, sources.use =c('FPR_Mo'), 
                 targets.use = c(1:10),  comparison = c(1, 2), angle.x = 45)






                       
