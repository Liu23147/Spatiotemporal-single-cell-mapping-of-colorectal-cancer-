NAT <- Load10X_Spatial(data.dir = localdir)
NAT <- NormalizeData(NAT)
NAT <- FindVariableFeatures(NAT)
NAT <- ScaleData(NAT)


SpatialFeaturePlot(NAT, features = c('CXCL10','CXCR3'),image.alpha=0 ,
                   pt.size.factor = 20)
SpatialFeaturePlot(P1, features = c('ANXA1'),image.alpha=0)
SpatialFeaturePlot(P1, features = c("RGS13"))

SpatialDimPlot(NAT, label = T, repel = T, label.size = 4)

NAT<-pre

NAT <- SketchData(
  object = NAT,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(NAT) <- "sketch"
NAT <- ScaleData(NAT)
NAT <- RunPCA(NAT, assay = "sketch", reduction.name = "pca.NAT.sketch", verbose = T)
NAT <- FindNeighbors(NAT, reduction = "pca.NAT.sketch", dims = 1:50)
NAT <- RunUMAP(NAT, reduction = "pca.NAT.sketch", reduction.name = "umap.NAT.sketch", return.model = T, dims = 1:50, verbose = T)

ref=readRDS(ref.RDS)
Idents(ref) <- "group"
ref<-subset(ref, idents=c('Pre'))

Idents(ref) <- "celltype"
counts <- ref[["SCT"]]$counts
cluster <- as.factor(ref$celltype)
nUMI <- as.integer(ref$nCount_RNA)
names(nUMI) <- colnames(counts)  # Add cell barcodes as names

levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# Ensure cluster also has cell barcode names (if not already)
names(cluster) <- colnames(counts)  # Only if cluster is a vector without names

reference <- Reference(counts, cluster, nUMI) 

counts_hd <- NAT[["sketch"]]$counts
NAT_cells_hd <- colnames(NAT[["sketch"]])
coords <- GetTissueCoordinates(NAT)[NAT_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

RCTD <- create.RCTD(query, reference, max_cores = 1,
                    CELL_MIN_INSTANCE = 0)

RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
NAT <- AddMetaData(NAT, metadata = RCTD@results$results_df)

NAT$first_type <- as.character(NAT$first_type)
NAT$first_type[is.na(NAT$first_type)] <- "Unknown"

NAT <- ProjectData(
  object = NAT,
  assay = "Spatial",
  full.reduction = "pca.NAT",
  sketched.assay = "sketch",
  sketched.reduction = "pca.NAT.sketch",
  umap.model = "umap.NAT.sketch",
  dims = 1:50,
  refdata = list(full_first_type= "first_type")
)

DefaultAssay(NAT) <- "Spatial"

# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
NAT[[]][, "full_first_type"] <- "Unknown"
NAT$full_first_type[Cells(NAT)] <- NAT$first_type[Cells(NAT)]

library(semla)
library(DT)
SANT <- UpdateSeuratForSemla(ANT)
library(DT)
res_net <-RunNeighborhoodEnrichmentTest(
  SANT,
  column_name='full_first_type',
  column_labels = c( "Fibroblast"     ,  "Epithelial cell" , "Endothelial cell" ,"IGA-P"    ,        "GPR34-Ma"    ,     "cDC3"     ,        "IGHD-B" ,         
                     "CD4-Naive-T" ,     "CD24-B" ,          "CD27-B"   ,        "SELENOP-Ma"   ,    "FPR-Mo"    ,       "CXCL13-T"  ,       "CD16-NK"  ,        "cDC1"  ,          
                     "CCL20-T"   ,       "IGM-P"   ,         "CXCL10-Ma"  ,      "cDC2" ,            "ANXA1-T" ,         "IGG-P"  ,          "SPP1-Ma"  ,        "IL17A-T"  ,       
                     "CD8-IFNG-T" ,      "LAG3-Tex"   ,      "RGS13-B"   ,       "FOXP3-Treg" ,      "LRMDA-Mo"    ,     "Neutrophil"    ,   "CCL3L3-Mo"   ,     "MT1H-Ma"  ,       
                     "NKT"   ,           "HSP-B"  ,          "Naive-T" ,         "MKI67-Ma"   ,      "HSPA1A-T"     ,    "GZMK-CD8-T"    ,   "pDC" ),
  n_permutations = 1000,
  nCores = parallel::detectCores() - 1,
  seed = 123,
  verbose = TRUE
)

datatable(
  res_net |> arrange(desc(abs(z_score))), 
  rownames = F, 
  caption = "Neighborhood enrichment analysis output"
)
cluster_labels <- paste0("Label_", sort(unique(SANT$full_first_type)))
hm_plot_data <- res_net

hm_plot_data$label_1 <- factor(hm_plot_data$label_1, levels = cluster_labels)
hm_plot_data$label_2 <- factor(hm_plot_data$label_2, levels = cluster_labels)

ggplot(hm_plot_data, aes(label_1, label_2, fill= z_score)) + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradient2(low = "#0474BA",
                       mid = "grey90",
                       high = "#F17720") +
  labs(x="", y="", title="Neighborhood enrichment", fill="Z-score") +
  theme_bw() +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust=0.5, size=12, face = "bold"),
        axis.text = element_text(size=10),
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10))

pre <- RunBanksy(pre,
                 lambda = 0.8, verbose = TRUE,
                 assay = "Spatial", slot = "data", features = "variable",
                 k_geom = 50
)
DefaultAssay(pre) <- "BANKSY"
pre <- RunPCA(pre, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(pre), npcs = 30)
pre <- FindNeighbors(pre, reduction = "pca.banksy", dims = 1:30)
pre <- FindClusters(pre, cluster.name = "banksy_cluster", resolution = 0.5)
Idents(pre) <- "banksy_cluster"
SpatialDimPlot(pre, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)

cells <- CellsByIdentities(pre, idents = c(18))

SpatialDimPlot(
  pre,
  pt.size.factor = 5,
  cells.highlight = cells,  
  cols.highlight = c('#4991C1',
                     'grey50'                 ),  
  facet.highlight = FALSE,  
  combine = TRUE)


SpatialDimPlot(pre, group.by = "banksy_cluster", label = F, repel = T, label.size = 4)+ SpatialDimPlot(P1, group.by = "banksy_cluster", label = F, repel = T, label.size = 4)

cells <- CellsByIdentities(P1, idents = c(18))

SpatialDimPlot(
  P1,
  pt.size.factor = 5,
  cells.highlight = cells, 
  cols.highlight = c('#C6307C',
                     'grey50'), 
  facet.highlight = FALSE,  
  combine = TRUE)
cluster1.markers <- FindMarkers(P1, ident.1 = '7', verbose = TRUE, logfc.threshold=0.05)













