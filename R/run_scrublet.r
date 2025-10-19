# Run doublet discrimination using Scrublet via reticulate
#' @export
run_scrublet <- function(seurat_obj, 
                                expected_doublet_rate = 0.06,
                                min_counts = 2,
                                min_cells = 3,
                                min_gene_variability_pctl = 85,
                                filter_doublets = TRUE,
                                verbose = TRUE) {
  
  # Load required library
  if (!require(reticulate, quietly = TRUE)) {
    stop("reticulate package is required. Please install it.")
  }
  
  # Import scrublet
  scrublet <- import("scrublet")
  
  if (verbose) {
    cat("Running Scrublet on", ncol(seurat_obj), "cells...\n")
  }
  
  # Get raw counts matrix (transpose to cells x genes for Scrublet)
  counts_matrix <- t(as.matrix(GetAssayData(seurat_obj, slot = "counts")))
  
  # Initialize Scrublet
  scrub <- scrublet$Scrublet(counts_matrix, expected_doublet_rate = expected_doublet_rate)
  
  # Calculate doublet scores
  result <- scrub$scrub_doublets(min_counts = min_counts, 
                                 min_cells = min_cells,
                                 min_gene_variability_pctl = min_gene_variability_pctl)
  
  doublet_scores <- result[[1]]
  predicted_doublets <- result[[2]]
  
  # Add results to Seurat object metadata
  seurat_obj$scrublet_score <- doublet_scores
  seurat_obj$scrublet_prediction <- predicted_doublets
  
  if (verbose) {
    n_doublets <- sum(predicted_doublets)
    doublet_rate <- round(n_doublets / length(predicted_doublets) * 100, 2)
    cat("  Detected", n_doublets, "doublets (", doublet_rate, "%)\n")
  }
  
  # Filter doublets if requested
  if (filter_doublets) {
    seurat_obj <- subset(seurat_obj, scrublet_prediction == FALSE)
    if (verbose) {
      cat("  Filtered object contains", ncol(seurat_obj), "singlets\n")
    }
  }
  
  return(seurat_obj)
}