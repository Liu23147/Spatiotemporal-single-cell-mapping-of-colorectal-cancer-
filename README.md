scVirtualKO: In Silico Knockout Analysis for Single-Cell RNA-seq

scVirtualKO is an R-based tool designed to perform in silico perturbation analysis on single-cell RNA-sequencing (scRNA-seq) data. By leveraging Random Forest regression at the single-cell level, it models the dependency between a target gene (KO gene) and downstream targets, allowing you to predict the shift in gene expression upon virtual knockout.

Key Idea: Train on the population, predict on single cells, and validate with multi-modal prior knowledge (PPI, TF motifs, Ligand-Receptor networks).

Cell-Level Modeling: Builds independent Random Forest models for each target gene using the KO gene and background context genes as predictors.

Virtual Perturbation: Simulates the effect of gene knockout (KO_Gene = 0) while maintaining the cellular context (background genes).

Multi-Modal Validation: automatically integrates prior knowledge to tag results:

TF: Transcription Factor motifs (via DoRothEA).

PPI: Protein-Protein Interactions (via OmniPath).

Ligand/Receptor: Downstream signaling targets (via NicheNet).

Smart Filtering: Automatically excludes mitochondrial (MT-) and ribosomal (RPS/RPL) genes to reduce noise.

Flexible Background Control: Supports automatic High-Variable Genes (HVGs) selection, manual background gene lists, or pure univariate modeling (no background).

#Ensure you have the following R packages installed:
install.packages(c("Seurat", "ranger", "dplyr", "tibble", "Matrix", "ggplot2", "ggrepel"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("OmnipathR", "dorothea", "nichenetr"))
library(Seurat)
source("scVirtualKO_function.R") # Load the function script

# 1. Load your Seurat object
seurat_obj <- readRDS("your_data.rds")
ppi_source=read.csv("ppi_source.csv")
niche_net=readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
niche_mat=readRDS(url("https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
# 2. Run Virtual KO
To accurately model the effect of a KO gene, it is crucial to provide a stable "background context." This helps the model distinguish between specific gene regulation and general cell states (e.g., cell cycle phase, metabolic activity).

my_background=c(GAPDH","ACTB","MKI67","TOP2A","PCNA","CDK1","CCNB1")

results <- scVirtualKO_Loop (scRNA_metacell , ko_gene = "FPR1", mode = c('Receptor'),manual_background_genes = my_background  ,
                            n_targets = 1000, 
                            n_trees = 2000,
                            seed = 123,
                            verbose = TRUE, 
                            nichenet_mat = niche_mat,  
                            ppi_source = ppi_source ,
                            lr_net = niche_net,
                            candidate_genes = NULL)

# 3. View Results
head(results)




#Parameters
seurat_obj	A Seurat object containing single-cell expression data.	Required

ko_gene	The gene symbol to virtually knockout.	Required

manual_background_genes	Control background context:- NULL: Auto-select top 50 HVGs.- FALSE: Disable background.- Vector: Use specific genes (Recommended).	NULL

confounders	A vector of metadata columns (e.g., "seurat_clusters", "pseudotime") to include as covariates. Essential for fixing cell identity during prediction.	NULL

mode	Validation strategies. Options: "TF", "PPI", "Ligand", "Receptor".	All

n_targets	Number of additional high-variance genes to predict.	100

candidate_genes	A vector of specific target genes you want to predict (optional).	NULL

ppi_source	Custom dataframe for PPI validation (columns: source_genesymbol, target_genesymbol).	NULL
