####nichenetr####
library(pheatmap)
library(nichenetr)
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj =scRNA10, 
  receiver = c('CXCL13-Th',  'FCGR3A-NK',  'GZMK-Teff',   'HSPA1A-T',     'IFNG-T',    'IL17-Th',    'Naive-T','Treg',    'XCL1-NK'), 
  condition_colname = "tissue", condition_oi = 'HPV-CC', condition_reference = 'NHPV-CC',geneset='up',assay_oi='SCT',
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
