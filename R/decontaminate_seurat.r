# =============================================================================
# DECONTAMINATE SEURAT FUNCTION
# =============================================================================

#' Remove ambient RNA contamination from Seurat objects using decontX
#' 
#' This function applies decontX decontamination to remove ambient RNA contamination
#' from single-cell RNA-seq data. Works with single Seurat objects or lists of objects.
#' 
#' @param seurat_obj Seurat object or named list of Seurat objects to decontaminate
#' @param assay Name of assay to decontaminate (default: "RNA")
#' @param background Background matrix specification. Can be:
#'   - NULL: decontX will estimate background automatically
#'   - Character: name of assay to use as background (e.g., "RAW")
#'   - Matrix: dgCMatrix or matrix to use directly as background
#' @param max_iter Maximum iterations for decontX (default: 500)
#' @param convergence Convergence threshold for decontX (default: 0.001)
#' @param save_plots Logical, whether to save diagnostic plots (default: TRUE)
#' @param plot_dir Directory to save plots (default: "decontX_plots/")
#' @param plot_prefix Optional prefix for plot filenames
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param BPPARAM BiocParallel parameter for parallel processing (optional, for lists only)
#' @param ... Additional arguments passed to decontX()
#' 
#' @return Single decontaminated Seurat object (if input was single object) or 
#'   named list of decontaminated Seurat objects (if input was list) with:
#'   - New "decontX" assay with decontaminated counts
#'   - Original assay preserved unchanged
#'   - Added metadata: contamination_score, decontX_clusters, etc.
#'   - Plots saved to specified directory
#' @export
decontaminate_seurat <- function(seurat_obj,
                                 assay = "RNA",
                                 background = NULL,
                                 max_iter = 500,
                                 convergence = 0.001,
                                 save_plots = TRUE,
                                 plot_dir = "decontX_plots/",
                                 plot_prefix = NULL,
                                 verbose = TRUE,
                                 BPPARAM = NULL,
                                 ...) {
  
  # Load required packages
  required_packages <- c("celda", "SingleCellExperiment", "ggplot2", "patchwork", "scales", "Seurat")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not available"))
    }
  }
  
  # Check if input is a list of Seurat objects or single object
  if (is.list(seurat_obj) && all(sapply(seurat_obj, function(x) inherits(x, "Seurat")))) {
    return(decontaminate_seurat_list(seurat_obj, assay, background, max_iter, convergence,
                                     save_plots, plot_dir, plot_prefix, verbose, 
                                     BPPARAM, ...))
  } else if (inherits(seurat_obj, "Seurat")) {
    return(decontaminate_seurat_single(seurat_obj, assay, background, max_iter, convergence,
                                       save_plots, plot_dir, plot_prefix, verbose, ...))
  } else {
    stop("Input must be a Seurat object or a list of Seurat objects")
  }
}

# =============================================================================
# SINGLE OBJECT PROCESSING FUNCTION
# =============================================================================

decontaminate_seurat_single <- function(seurat_obj,
                                        assay = "RNA",
                                        background = NULL,
                                        max_iter = 500,
                                        convergence = 0.001,
                                        save_plots = TRUE,
                                        plot_dir = "decontX_plots/",
                                        plot_prefix = NULL,
                                        verbose = TRUE,
                                        ...) {
  
  if (verbose) {
    message("Starting decontX decontamination...")
  }
  
  # Get sample name
  sample_name <- seurat_obj@project.name
  
  # Get the count matrix to decontaminate
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste("Assay", assay, "not found in Seurat object"))
  }
  
  counts_data <- GetAssayData(seurat_obj, assay = assay, slot = "counts")
  
  if (verbose) {
    message(paste("Processing sample:", sample_name))
    message(paste("Decontaminating assay:", assay))
    message(paste("Input matrix:", nrow(counts_data), "genes x", ncol(counts_data), "cells"))
  }
  
  # Prepare background matrix
  background_matrix <- prepare_background_matrix(seurat_obj, assay, background, verbose)
  
  # Run decontX
  if (verbose) message("Running decontX...")
  
  decontX_results <- run_decontx(counts_data, background_matrix, max_iter, convergence, verbose, ...)
  
  # Extract results
  decontaminated_counts <- decontX_results$decontaminated_counts
  contamination_scores <- decontX_results$contamination_scores
  decontx_clusters <- decontX_results$clusters
  umap_coords <- decontX_results$umap_coords
  
  if (verbose) {
    message(paste("Original matrix:", nrow(counts_data), "genes x", ncol(counts_data), "cells"))
    message(paste("Decontaminated matrix:", nrow(decontaminated_counts), "genes x", ncol(decontaminated_counts), "cells"))
    
    # Debug information
    message("Checking matrix compatibility...")
    message(paste("Original cell names: first 3 =", paste(head(colnames(counts_data), 3), collapse = ", ")))
    message(paste("Decontaminated cell names: first 3 =", paste(head(colnames(decontaminated_counts), 3), collapse = ", ")))
    message(paste("Original gene names: first 3 =", paste(head(rownames(counts_data), 3), collapse = ", ")))
    message(paste("Decontaminated gene names: first 3 =", paste(head(rownames(decontaminated_counts), 3), collapse = ", ")))
  }
  
  # Ensure the decontaminated matrix has proper dimnames
  if (is.null(colnames(decontaminated_counts))) {
    colnames(decontaminated_counts) <- colnames(counts_data)
    if (verbose) message("Added missing column names to decontaminated matrix")
  }
  
  if (is.null(rownames(decontaminated_counts))) {
    # If genes were filtered, we need to match the remaining genes
    if (nrow(decontaminated_counts) == nrow(counts_data)) {
      rownames(decontaminated_counts) <- rownames(counts_data)
    } else {
      # DecontX may have filtered some genes, create generic names
      rownames(decontaminated_counts) <- paste0("Gene_", 1:nrow(decontaminated_counts))
      if (verbose) message("Warning: Gene filtering occurred, using generic gene names")
    }
    if (verbose) message("Added missing row names to decontaminated matrix")
  }
  
  # Ensure matrix is in the right format (sparse matrix)
  if (!inherits(decontaminated_counts, "sparseMatrix")) {
    if (verbose) message("Converting to sparse matrix format...")
    decontaminated_counts <- Matrix::Matrix(decontaminated_counts, sparse = TRUE)
  }
  
  # Create new assay with decontaminated counts instead of replacing original
  if (verbose) message("Creating decontX assay...")
  
  tryCatch({
    seurat_obj[["decontX"]] <- CreateAssay5Object(counts = decontaminated_counts)
    if (verbose) message("Successfully created decontX assay")
  }, error = function(e) {
    message("Error creating decontX assay:", e$message)
    message("Attempting to fix matrix format...")
    
    # Try to fix common issues
    # Ensure dimensions match expectations
    if (ncol(decontaminated_counts) != ncol(counts_data)) {
      stop("Cell count mismatch between original and decontaminated matrices")
    }
    
    # Try creating with explicit column names matching the seurat object
    current_cells <- colnames(seurat_obj)
    if (length(current_cells) == ncol(decontaminated_counts)) {
      colnames(decontaminated_counts) <- current_cells
      if (verbose) message("Matched cell names to current Seurat object")
    }
    
    # Try again
    seurat_obj[["decontX"]] <<- CreateAssay5Object(counts = decontaminated_counts)
    if (verbose) message("Successfully created decontX assay after fixing")
  })
  
  # Set decontX as default assay
  DefaultAssay(seurat_obj) <- "decontX"
  
  # Add metadata
  seurat_obj$contamination_score <- contamination_scores
  seurat_obj$decontX_clusters <- decontx_clusters
  
  # Add UMAP coordinates if available
  if (!is.null(umap_coords) && is.matrix(umap_coords) && ncol(umap_coords) >= 2) {
    seurat_obj$decontX_UMAP_1 <- umap_coords[,1]
    seurat_obj$decontX_UMAP_2 <- umap_coords[,2]
    if (verbose) message("Added decontX UMAP coordinates to metadata")
  }
  
  # Print summary
  if (verbose) {
    print_decontx_summary(contamination_scores, sample_name)
  }
  
  # Create plots
  if (save_plots) {
    create_decontx_plots_enhanced(
      original_counts = counts_data,
      decontaminated_counts = decontaminated_counts,
      contamination_scores = contamination_scores,
      umap_coords = umap_coords,
      sample_name = sample_name,
      plot_dir = plot_dir,
      plot_prefix = plot_prefix,
      verbose = verbose
    )
  }
  
  if (verbose) {
    message("DecontX decontamination completed successfully!")
    message("Original assay preserved, decontaminated counts stored in 'decontX' assay")
    message("Default assay set to 'decontX'")
  }
  
  return(seurat_obj)
}

# =============================================================================
# LIST PROCESSING FUNCTION
# =============================================================================

decontaminate_seurat_list <- function(seurat_list,
                                      assay = "RNA",
                                      background = NULL,
                                      max_iter = 500,
                                      convergence = 0.001,
                                      save_plots = TRUE,
                                      plot_dir = "decontX_plots/",
                                      plot_prefix = NULL,
                                      verbose = TRUE,
                                      BPPARAM = NULL,
                                      ...) {
  
  if (verbose) {
    message(paste("Processing", length(seurat_list), "Seurat objects for decontamination"))
  }
  
  # Get sample names
  sample_names <- names(seurat_list)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_along(seurat_list))
    names(seurat_list) <- sample_names
    if (verbose) {
      message("No names found, using default names: ", paste(sample_names, collapse = ", "))
    }
  }
  
  # Internal processing function for one object
  process_one_object <- function(i) {
    obj <- seurat_list[[i]]
    name <- sample_names[i]
    
    if (verbose) {
      message(paste("Processing sample:", name))
    }
    
    # Create sample-specific plot prefix if provided
    current_plot_prefix <- if (!is.null(plot_prefix)) {
      paste0(plot_prefix, name, "_")
    } else {
      paste0(name, "_")
    }
    
    # Process the object
    decontaminated_obj <- decontaminate_seurat_single(
      seurat_obj = obj,
      assay = assay,
      background = background,
      max_iter = max_iter,
      convergence = convergence,
      save_plots = save_plots,
      plot_dir = plot_dir,
      plot_prefix = current_plot_prefix,
      verbose = verbose,
      ...
    )
    
    return(decontaminated_obj)
  }
  
  # Process objects (parallel or serial)
  if (!is.null(BPPARAM)) {
    # Parallel processing
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      warning("BiocParallel not available, falling back to serial processing")
      decontaminated_objects <- lapply(seq_along(seurat_list), process_one_object)
    } else {
      if (verbose) {
        message("Processing objects in parallel...")
      }
      
      decontaminated_objects <- tryCatch({
        BiocParallel::bplapply(seq_along(seurat_list), process_one_object, BPPARAM = BPPARAM)
      }, error = function(e) {
        warning(paste("Parallel processing failed:", e$message, 
                      "\nFalling back to serial processing"))
        lapply(seq_along(seurat_list), process_one_object)
      })
    }
  } else {
    # Serial processing
    if (verbose) {
      message("Processing objects serially...")
    }
    decontaminated_objects <- lapply(seq_along(seurat_list), process_one_object)
  }
  
  # Restore names
  names(decontaminated_objects) <- sample_names
  
  if (verbose) {
    message(paste("Successfully decontaminated", length(decontaminated_objects), "Seurat objects"))
    
    # Print summary for all objects
    cat("\n=== Decontamination Summary for All Samples ===\n")
    for (i in seq_along(decontaminated_objects)) {
      obj <- decontaminated_objects[[i]]
      name <- names(decontaminated_objects)[i]
      mean_contamination <- round(mean(obj$contamination_score, na.rm = TRUE), 4)
      median_contamination <- round(median(obj$contamination_score, na.rm = TRUE), 4)
      high_contamination <- sum(obj$contamination_score > 0.1, na.rm = TRUE)
      
      cat(sprintf("%s: Mean contamination: %.4f, Median: %.4f, High (>10%%): %d cells\n", 
                  name, mean_contamination, median_contamination, high_contamination))
    }
    cat("\n")
  }
  
  return(decontaminated_objects)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to prepare background matrix based on different input types
prepare_background_matrix <- function(seurat_obj, assay, background, verbose) {
  
  if (is.null(background)) {
    if (verbose) message("No background specified - decontX will estimate background automatically")
    return(NULL)
  }
  
  if (is.character(background)) {
    # Background is an assay name
    if (!background %in% names(seurat_obj@assays)) {
      stop(paste("Background assay", background, "not found in Seurat object"))
    }
    
    # Get barcodes from both assays
    active_barcodes <- colnames(GetAssayData(seurat_obj, assay = assay))
    background_barcodes <- colnames(GetAssayData(seurat_obj, assay = background))
    
    # Find cells in background but not in active assay
    empty_droplet_barcodes <- setdiff(background_barcodes, active_barcodes)
    
    if (length(empty_droplet_barcodes) == 0) {
      warning("No empty droplets found - all cells in background assay are also in active assay")
      return(NULL)
    }
    
    background_matrix <- GetAssayData(seurat_obj, assay = background)[, empty_droplet_barcodes, drop = FALSE]
    
    if (verbose) {
      message(paste("Using", ncol(background_matrix), "empty droplets from assay", background, "as background"))
    }
    
    return(background_matrix)
    
  } else if (is.matrix(background) || inherits(background, "sparseMatrix")) {
    # Background is a matrix
    if (verbose) {
      message(paste("Using provided matrix with", ncol(background), "cells as background"))
    }
    return(background)
    
  } else {
    stop("Background must be NULL, an assay name (character), or a count matrix")
  }
}

# Function to run decontX and extract results
run_decontx <- function(counts_data, background_matrix, max_iter, convergence, verbose, ...) {
  
  # Run decontX
  if (!is.null(background_matrix)) {
    decontX_results <- decontX(
      x = counts_data,
      background = background_matrix,
      maxIter = max_iter,
      convergence = convergence,
      verbose = verbose,
      ...
    )
  } else {
    decontX_results <- decontX(
      x = counts_data,
      maxIter = max_iter,
      convergence = convergence,
      verbose = verbose,
      ...
    )
  }
  
  # Extract results based on return type
  if (is(decontX_results, "SingleCellExperiment")) {
    decontaminated_counts <- assay(decontX_results, "decontXcounts")
    contamination_scores <- colData(decontX_results)$decontX_contamination
    clusters <- colData(decontX_results)$decontX_clusters
  } else if (is.list(decontX_results)) {
    decontaminated_counts <- decontX_results$decontXcounts
    contamination_scores <- decontX_results$contamination
    clusters <- decontX_results$z
  } else {
    stop("Unexpected decontX results format")
  }
  
  # Extract UMAP coordinates if available
  umap_coords <- NULL
  if (is.list(decontX_results) && 
      !is.null(decontX_results$estimates) && 
      !is.null(decontX_results$estimates$all_cells) && 
      !is.null(decontX_results$estimates$all_cells$UMAP)) {
    umap_coords <- decontX_results$estimates$all_cells$UMAP
  }
  
  return(list(
    decontaminated_counts = decontaminated_counts,
    contamination_scores = contamination_scores,
    clusters = clusters,
    umap_coords = umap_coords
  ))
}

# Function to print decontX summary
print_decontx_summary <- function(contamination_scores, sample_name) {
  cat("\n=== DecontX Results Summary ===\n")
  cat("Sample:", sample_name, "\n")
  cat("Mean contamination score:", round(mean(contamination_scores, na.rm = TRUE), 4), "\n")
  cat("Median contamination score:", round(median(contamination_scores, na.rm = TRUE), 4), "\n")
  cat("Max contamination score:", round(max(contamination_scores, na.rm = TRUE), 4), "\n")
  cat("Cells with >10% contamination:", sum(contamination_scores > 0.1, na.rm = TRUE), "\n")
  cat("Cells with >25% contamination:", sum(contamination_scores > 0.25, na.rm = TRUE), "\n\n")
}

# Enhanced function to create decontX diagnostic plots
create_decontx_plots_enhanced <- function(original_counts, decontaminated_counts, contamination_scores,
                                          umap_coords, sample_name, plot_dir, plot_prefix, verbose) {
  
  # Create plots directory if needed
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Plot 1: Contamination score distribution
  df_contamination <- data.frame(
    contamination = contamination_scores,
    cell_id = seq_along(contamination_scores)
  )
  
  p1 <- ggplot(df_contamination, aes(x = contamination)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = mean(contamination_scores, na.rm = TRUE), color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = median(contamination_scores, na.rm = TRUE), color = "blue", linetype = "dashed", size = 1) +
    labs(
      x = "Contamination Score",
      y = "Number of Cells",
      title = paste("Contamination Score Distribution -", sample_name),
      subtitle = paste("Mean:", round(mean(contamination_scores, na.rm = TRUE), 3), 
                       "| Median:", round(median(contamination_scores, na.rm = TRUE), 3))
    ) +
    theme_minimal()
  
  # Plot 2: UMAP with contamination scores (if available)
  p2 <- NULL
  if (!is.null(umap_coords) && is.matrix(umap_coords) && ncol(umap_coords) >= 2) {
    umap_df <- data.frame(
      UMAP_1 = umap_coords[,1],
      UMAP_2 = umap_coords[,2],
      contamination = contamination_scores
    )
    
    p2 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = contamination)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_gradient(low = "blue", high = "red", name = "Contamination\nScore") +
      labs(
        x = "UMAP 1",
        y = "UMAP 2",
        title = paste("UMAP with Contamination Scores -", sample_name)
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  }
  
  # Plot 3: Contamination vs total UMI count
  total_umis <- colSums(original_counts)
  df_umi <- data.frame(
    total_UMI = total_umis,
    contamination = contamination_scores
  )
  
  p3 <- ggplot(df_umi, aes(x = total_UMI, y = contamination)) +
    geom_point(alpha = 0.6, size = 0.5) +
    geom_smooth(method = "loess", color = "red", se = TRUE) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(
      x = "Total UMI Count",
      y = "Contamination Score",
      title = paste("Contamination vs UMI Count -", sample_name)
    ) +
    theme_minimal()
  
  # Plot 4: Gene expression changes (before vs after)
  gene_vars_before <- apply(original_counts, 1, var)
  top_var_genes <- names(sort(gene_vars_before, decreasing = TRUE))[1:min(1000, length(gene_vars_before))]
  
  total_before <- rowSums(original_counts[top_var_genes, ])
  total_after <- rowSums(decontaminated_counts[top_var_genes, ])
  
  df_comparison <- data.frame(
    gene = names(total_before),
    before = total_before,
    after = total_after
  )
  
  p4 <- ggplot(df_comparison, aes(x = before, y = after)) +
    geom_point(alpha = 0.6, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(
      x = "Total Expression Before DecontX",
      y = "Total Expression After DecontX",
      title = paste("Gene Expression Changes -", sample_name)
    ) +
    theme_minimal()
  
  # Create combined plot layout
  if (!is.null(p2)) {
    combined_plot <- (p1 | p2) / (p3 | p4)
    fig_height <- 8
  } else {
    combined_plot <- p1 / (p3 | p4)
    fig_height <- 6
  }
  
  # Save plot
  if(is.null(plot_prefix)) {
    plot_filename <- file.path(plot_dir, paste0(sample_name, "_decontX_analysis.tiff"))
  } else {
    plot_filename <- file.path(plot_dir, paste0(plot_prefix, sample_name, "_decontX_analysis.tiff"))
  }
  
  ggsave(plot_filename, combined_plot, width = 12, height = fig_height, dpi = 300)
  
  if (verbose) {
    message("DecontX analysis plots saved to:", plot_filename)
  }
  
  return(plot_filename)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Single object with automatic background estimation
# decontaminated_seurat <- decontaminate_seurat(seurat_obj)

# Single object with background from RAW assay (typical after filter_empty_cells)  
# decontaminated_seurat <- decontaminate_seurat(seurat_obj, background = "RAW")

# Single object with custom background matrix
# background_matrix <- some_empty_droplet_matrix
# decontaminated_seurat <- decontaminate_seurat(seurat_obj, background = background_matrix)

# Access different count matrices after decontamination:
# original_counts <- GetAssayData(decontaminated_seurat, assay = "RNA", slot = "counts")     # Original filtered
# decontaminated_counts <- GetAssayData(decontaminated_seurat, assay = "decontX", slot = "counts")  # Decontaminated
# raw_counts <- GetAssayData(decontaminated_seurat, assay = "RAW", slot = "counts")         # Raw unfiltered (if available)

# Check contamination scores
# summary(decontaminated_seurat$contamination_score)
# table(decontaminated_seurat$decontX_clusters)

# List of objects (serial)
# seurat_list <- list(Sample1 = obj1, Sample2 = obj2, Sample3 = obj3)
# decontaminated_list <- decontaminate_seurat(seurat_list, background = "RAW")

# List of objects (parallel)
# library(BiocParallel)
# decontaminated_list <- decontaminate_seurat(
#   seurat_list, 
#   background = "RAW",
#   BPPARAM = MulticoreParam(workers = 4)
# )

# Custom parameters
# decontaminated_list <- decontaminate_seurat(
#   seurat_list,
#   assay = "RNA",
#   background = "RAW",
#   max_iter = 1000,
#   convergence = 0.0001,
#   plot_dir = "my_decontX_plots/",
#   plot_prefix = "experiment1_",
#   BPPARAM = MulticoreParam(workers = 2)
# )