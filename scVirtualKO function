library(Seurat)
library(ranger)
library(dplyr)
library(tibble)
library(Matrix)

# ==============================================================================
# 1. Main Function: scVirtualKO
# ==============================================================================

#' @title In Silico Knockout Analysis (Looping/Univariate Approach)
#' @description Models gene dependencies per cell to predict knockout effects.
#' @param seurat_obj Seurat Object.
#' @param ko_gene Gene symbol to knockout.
#' @param candidate_genes (Optional) Manual list of TARGET genes to analyze.
#' @param manual_background_genes (Optional) 
#'        - NULL (Default): Auto-select Top 50 HVGs as background.
#'        - FALSE: DISABLE background genes (Model: Target ~ KO_Gene only).
#'        - Character Vector: Use specific genes as background.
#' @param use_all_genes (Bool) If TRUE, selects automatic background from all genes.
#' @param mode Validation modes: "TF", "Ligand", "Receptor", "PPI".
#' @param nichenet_mat NicheNet Ligand-Target Matrix.
#' @param lr_net NicheNet Ligand-Receptor Network.
#' @param ppi_source Manual dataframe for PPI.
#' @param n_targets Number of automatic targets to add.
#' @param n_trees Number of trees for Random Forest.
#' @param seed Random seed.
#' @param verbose Print progress messages.
scVirtualKO <- function(seurat_obj, 
                             ko_gene, 
                             candidate_genes = NULL,
                             manual_background_genes = NULL, 
                             use_all_genes = FALSE,
                             mode = c("Raw", "TF", "Ligand", "Receptor", "PPI"), 
                             nichenet_mat = NULL,  
                             lr_net = NULL,
                             ppi_source = NULL,    
                             n_targets = 100, 
                             n_trees = 500,
                             seed = 123,
                             verbose = TRUE) {
  
  mode <- match.arg(mode, several.ok = TRUE)
  set.seed(seed)
  
  if (verbose) message(paste0(">>> [Main] Analyzing KO Gene: ", ko_gene))
  if (verbose) message(paste0(">>> [Mode] Active validation modes: ", paste(mode, collapse = ", ")))
  
  if (!ko_gene %in% rownames(seurat_obj)) {
    stop(paste0("Error: Gene ", ko_gene, " not found in Seurat object."))
  }
  
  # --- Step 1: Run Core Looping Prediction ---
  if (identical(manual_background_genes, FALSE)) {
    bg_msg <- "DISABLED (Target ~ KO_Gene only)"
  } else if (!is.null(manual_background_genes)) {
    bg_msg <- "User-Defined Manual List"
  } else {
    bg_msg <- ifelse(use_all_genes, "All Genes (Auto)", "HVGs (Auto)")
  }
  
  if (verbose) message(paste0(">>> [Step 1] Starting loop modeling..."))
  if (verbose) message(paste0("    Filtering: Removing MT-, RPS, RPL genes from pool.")) # Log info
  if (verbose) message(paste0("    Background Source: ", bg_msg))
  
  rf_results <- .run_rf_loop_core(seurat_obj, ko_gene, candidate_genes, manual_background_genes, 
                                  use_all_genes, n_targets, n_trees, verbose)
  
  # --- Step 2: Sequential Validation ---
  rf_results$Validation_Tags <- ""
  # [Modified]: Removed Weight_Factor initialization
  
  # A. TF Validation
  if ("TF" %in% mode) {
    if (verbose) message(">>> [Validation] Checking TF Motifs (DoRothEA)...")
    tf_targets <- .get_tf_targets(ko_gene)
    hits <- rf_results$Gene %in% tf_targets
    rf_results$Validation_Tags[hits] <- paste(rf_results$Validation_Tags[hits], "Motif", sep = "; ")
    # [Modified]: Removed Weight calculation
  }
  
  # B. Ligand Validation
  if ("Ligand" %in% mode) {
    if (verbose) message(">>> [Validation] Checking Downstream Pathway (NicheNet Ligand)...")
    nn_targets <- .get_nichenet_targets(ko_gene, type="Ligand", nichenet_mat, lr_net=NULL)
    hits <- rf_results$Gene %in% nn_targets
    rf_results$Validation_Tags[hits] <- paste(rf_results$Validation_Tags[hits], "Pathway(L)", sep = "; ")
  }
  
  # C. Receptor Validation
  if ("Receptor" %in% mode) {
    if (verbose) message(">>> [Validation] Checking Downstream Pathway (NicheNet Receptor)...")
    if (is.null(lr_net)) stop("Error: 'lr_net' is required for Receptor mode.")
    nn_targets <- .get_nichenet_targets(ko_gene, type="Receptor", nichenet_mat, lr_net)
    hits <- rf_results$Gene %in% nn_targets
    rf_results$Validation_Tags[hits] <- paste(rf_results$Validation_Tags[hits], "Pathway(R)", sep = "; ")
  }
  
  # D. PPI Validation
  if ("PPI" %in% mode) {
    if (verbose) message(">>> [Validation] Checking Physical Interactions (PPI)...")
    ppi_targets <- .get_ppi_targets(ko_gene, ppi_source)
    hits <- rf_results$Gene %in% ppi_targets
    rf_results$Validation_Tags[hits] <- paste(rf_results$Validation_Tags[hits], "PPI", sep = "; ")
  }
  
  # --- Step 3: Formatting Output ---
  rf_results$Validation <- gsub("^; ", "", rf_results$Validation_Tags)
  rf_results$Validation <- trimws(rf_results$Validation)
  rf_results$Validation[rf_results$Validation == ""] <- "Predicted-Only"
  
  # [Modified]: Removed Hybrid_Score calculation
  
  if ("P_Value" %in% colnames(rf_results)) {
    rf_results$P_Adj <- p.adjust(rf_results$P_Value, method = "fdr")
  }
  
  if (!is.null(candidate_genes)) {
    rf_results$Source <- ifelse(rf_results$Gene %in% candidate_genes, "Candidate", "Auto")
  } else {
    rf_results$Source <- "Auto"
  }
  
  # [Modified]: Removed "Hybrid_Score" from cols_order
  cols_order <- c("Gene", "Source", "Delta_Mean", "Z_Score", "P_Value", "P_Adj", "KO_Importance", "Validation")
  cols_to_select <- intersect(cols_order, colnames(rf_results))
  extra_cols <- setdiff(colnames(rf_results), cols_to_select)
  
  # [Modified]: Also define extra cols to exclude intermediate tags
  extra_cols <- setdiff(extra_cols, c("Validation_Tags"))
  
  final_results <- rf_results[, c(cols_to_select, extra_cols)]
  
  # [Modified]: Sorting by absolute Z_Score (Magnitude of change) instead of Hybrid_Score
  final_results <- final_results[order(abs(final_results$Z_Score), decreasing = TRUE), ]
  rownames(final_results) <- NULL
  
  return(final_results)
}


# ==============================================================================
# 2. Helper Functions (Unchanged)
# ==============================================================================
.get_tf_targets <- function(ko_gene) {
  if (!requireNamespace("dorothea", quietly = TRUE)) return(character(0))
  data(dorothea_hs, package = "dorothea", envir = environment())
  dorothea_hs %>% filter(tf == ko_gene, confidence %in% c("A", "B", "C")) %>% pull(target)
}

.get_ppi_targets <- function(ko_gene, ppi_manual_data) {
  interactions <- NULL
  if (!is.null(ppi_manual_data)) {
    interactions <- ppi_manual_data
  } else if (requireNamespace("OmnipathR", quietly = TRUE)) {
    tryCatch({ interactions <- OmnipathR::import_all_interactions(organism = 9606) }, error = function(e) NULL)
  }
  if (is.null(interactions)) return(character(0))
  cols <- colnames(interactions)
  src_col <- grep("source.*symbol", cols, value = TRUE, ignore.case = TRUE)[1]
  tgt_col <- grep("target.*symbol", cols, value = TRUE, ignore.case = TRUE)[1]
  if (is.na(src_col) || is.na(tgt_col)) return(character(0))
  unique(c(interactions %>% filter(!!sym(src_col) == ko_gene) %>% pull(!!sym(tgt_col)),
           interactions %>% filter(!!sym(tgt_col) == ko_gene) %>% pull(!!sym(src_col))))
}

.get_nichenet_targets <- function(ko_gene, type, niche_mat, lr_net) {
  if (is.null(niche_mat)) return(character(0))
  query_ligands <- character(0)
  if (type == "Ligand") {
    if (ko_gene %in% rownames(niche_mat)) query_ligands <- ko_gene
  } else {
    upstream <- lr_net %>% filter(to == ko_gene) %>% pull(from) %>% unique()
    query_ligands <- intersect(upstream, rownames(niche_mat))
  }
  if (length(query_ligands) == 0) return(character(0))
  if (length(query_ligands) == 1) { scores <- niche_mat[query_ligands, ] } else { scores <- apply(niche_mat[query_ligands, , drop = FALSE], 2, max) }
  valid_scores <- scores[scores > 0]
  if(length(valid_scores) == 0) return(character(0))
  threshold <- quantile(valid_scores, 0.75) 
  names(valid_scores[valid_scores > threshold])
}

# ==============================================================================
# 3. Core Module: .run_rf_loop_core (With Filter)
# ==============================================================================

.run_rf_loop_core <- function(seurat_obj, ko_gene, candidate_genes, manual_background_genes, use_all_genes, n_targets, n_trees, verbose) {
  
  # --- A. Define Initial Pool ---
  if (use_all_genes) {
    pool <- rownames(seurat_obj)
  } else {
    pool <- VariableFeatures(seurat_obj)
    if (length(pool) < 50) {
      if(verbose) message("   Note: Insufficient variable features, recalculating Top 2000...")
      seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
      pool <- VariableFeatures(seurat_obj)
    }
  }
  
  # --- B. FILTER: Remove MT-, RPS, RPL Genes ---
  # Regex pattern for Mito and Ribosomal genes
  noise_pattern <- "^MT-|^RPS|^RPL"
  
  # Identify genes to remove (BUT keep the KO Gene itself if it happens to be one)
  genes_to_remove <- grep(noise_pattern, pool, value = TRUE)
  genes_to_remove <- setdiff(genes_to_remove, ko_gene) 
  
  if (length(genes_to_remove) > 0) {
    # Update the pool by removing these genes
    pool <- setdiff(pool, genes_to_remove)
    # if(verbose) message(paste0("   Removed ", length(genes_to_remove), " MT/RPS/RPL genes from modeling pool."))
  }
  
  # --- C. Handle Candidate Genes ---
  manual_targets <- character(0)
  if (!is.null(candidate_genes)) {
    valid_candidates <- intersect(candidate_genes, rownames(seurat_obj))
    
    # Filter candidate genes as well (optional, but consistent)
    cand_noise <- grep(noise_pattern, valid_candidates, value = TRUE)
    if(length(cand_noise) > 0) {
      if(verbose) message(paste0("   Warning: ", length(cand_noise), " candidate genes match MT/RPS/RPL pattern and will be excluded."))
      valid_candidates <- setdiff(valid_candidates, cand_noise)
    }
    
    manual_targets <- setdiff(valid_candidates, ko_gene)
    
    invalid <- setdiff(candidate_genes, rownames(seurat_obj))
    if (length(invalid) > 0) warning(paste("   Ignoring missing candidate genes:", paste(invalid, collapse=", ")))
  }
  
  # --- D. Define Final Targets ---
  remaining_pool <- setdiff(pool, c(ko_gene, manual_targets))
  n_auto <- min(n_targets, length(remaining_pool))
  auto_targets <- head(remaining_pool, n_auto)
  final_targets <- unique(c(manual_targets, auto_targets))
  
  # --- E. Define Background Predictors ---
  background_genes <- character(0)
  use_auto_bg <- FALSE
  
  if (identical(manual_background_genes, FALSE)) {
    # Disabled
  } else if (!is.null(manual_background_genes)) {
    valid_bg <- intersect(manual_background_genes, rownames(seurat_obj))
    
    # Also filter MT/RPS/RPL from manual background if user provided them
    bg_noise <- grep(noise_pattern, valid_bg, value = TRUE)
    if(length(bg_noise) > 0) {
      valid_bg <- setdiff(valid_bg, bg_noise)
      if(verbose) message(paste0("   Filtered ", length(bg_noise), " MT/RPS/RPL genes from manual background list."))
    }
    
    background_genes <- setdiff(valid_bg, c(ko_gene, final_targets))
    
    if (length(background_genes) == 0) {
      warning("   Warning: User-defined background list empty after filtering. Switching to AUTO mode.")
      use_auto_bg <- TRUE
    } else {
      if (verbose) message(paste0("   Using ", length(background_genes), " manual background genes."))
    }
  } else {
    use_auto_bg <- TRUE
  }
  
  if (use_auto_bg) {
    # Select from filtered pool
    bg_pool <- setdiff(pool, c(ko_gene, final_targets))
    n_bg <- min(50, length(bg_pool))
    background_genes <- if (n_bg > 0) head(bg_pool, n_bg) else character(0)
  }
  
  # --- F. Extract Data ---
  features_to_fetch <- unique(c(ko_gene, final_targets, background_genes))
  
  expr_data <- tryCatch({
    GetAssayData(seurat_obj, layer = "data")[features_to_fetch, ]
  }, error = function(e) {
    GetAssayData(seurat_obj, slot = "data")[features_to_fetch, ]
  })
  
  mat <- as.matrix(t(expr_data))
  
  # --- G. Clean Column Names ---
  original_colnames <- colnames(mat)
  clean_colnames <- make.names(original_colnames)
  name_map <- setNames(clean_colnames, original_colnames)
  colnames(mat) <- clean_colnames
  
  ko_safe <- name_map[[ko_gene]]
  if (is.null(ko_safe) || !ko_safe %in% colnames(mat)) stop(paste("Error: Cannot locate KO Gene column:", ko_gene))
  
  # --- H. Loop Prediction ---
  results_list <- list()
  if (verbose) message(paste0("   Training ", length(final_targets), " models..."))
  pb <- txtProgressBar(min = 0, max = length(final_targets), style = 3)
  
  for (i in seq_along(final_targets)) {
    target_orig <- final_targets[i]
    target_safe <- name_map[[target_orig]]
    
    if (is.na(target_safe)) next 
    
    bg_safe <- name_map[background_genes]
    bg_safe <- bg_safe[!is.na(bg_safe)]
    current_predictors <- unique(c(ko_safe, bg_safe))
    
    cols_to_use <- unique(c(target_safe, current_predictors))
    cols_to_use <- intersect(cols_to_use, colnames(mat))
    
    if (!target_safe %in% cols_to_use || !ko_safe %in% cols_to_use) next
    
    df_train <- as.data.frame(mat[, cols_to_use, drop = FALSE])
    
    rf_formula <- as.formula(paste(target_safe, "~ ."))
    
    fit <- tryCatch({
      ranger::ranger(formula = rf_formula, data = df_train, num.trees = n_trees, 
                     importance = "impurity", num.threads = 1, verbose = FALSE)
    }, error = function(e) return(NULL))
    
    if (is.null(fit)) next
    
    imp_val <- fit$variable.importance[ko_safe]
    if(is.na(imp_val)) imp_val <- 0
    
    df_ko <- df_train
    df_ko[, ko_safe] <- 0
    pred_wt <- predict(fit, data = df_train)$predictions
    pred_ko <- predict(fit, data = df_ko)$predictions
    
    delta <- pred_wt - pred_ko
    delta_mean <- mean(delta)
    
    p_val <- 1
    if (sd(delta) > 1e-10) {
      try({
        test_res <- t.test(pred_wt, pred_ko, paired = TRUE)
        p_val <- test_res$p.value
      }, silent = TRUE)
    }
    
    z_score <- delta_mean / (sd(delta) + 1e-6)
    
    results_list[[i]] <- data.frame(
      Gene = target_orig, 
      Delta_Mean = delta_mean, 
      Z_Score = z_score, 
      P_Value = p_val, 
      KO_Importance = imp_val, 
      stringsAsFactors = FALSE
    )
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  if (length(results_list) == 0) stop("Error: No results generated. (Check if filtering removed all targets)")
  
  res_df <- do.call(rbind, results_list)
  
  if (max(res_df$KO_Importance) > 0) {
    res_df$KO_Importance_Norm <- res_df$KO_Importance / max(res_df$KO_Importance)
  } else {
    res_df$KO_Importance_Norm <- 0
  }
  
  return(res_df)
}
