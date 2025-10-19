# =============================================================================
# FILTER EMPTY CELLS FUNCTION
# =============================================================================

#' Filter empty droplets from Seurat objects using emptyDrops
#' 
#' This function applies emptyDrops filtering to remove empty droplets from
#' single-cell RNA-seq data. It preserves background information for downstream
#' decontamination and adds relevant metadata to the Seurat object.
#' 
#' @param seurat_obj Seurat object to filter
#' @param threshold Minimum UMI count threshold for emptyDrops testing (default: 100)
#' @param test.ambient Logical, whether to test ambient contamination (default: TRUE)
#' @param FDR False discovery rate threshold for keeping cells (default: 0.001)
#' @param save_plots Logical, whether to save diagnostic plots (default: TRUE)
#' @param plot_dir Directory to save plots (default: "empty_cells_plots/")
#' @param plot_prefix Optional prefix for plot filenames
#' @param return_background Logical, whether to return background matrix for decontamination (default: TRUE)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param ... Additional arguments passed to emptyDrops()
#' 
#' @return Filtered Seurat object with:
#'   - RNA assay: Filtered count matrix (cells passing emptyDrops)
#'   - RAW assay: Original unfiltered count matrix (for background estimation)
#'   - Added metadata: cell_probability, barcode_rank, etc.
#'   - Plots saved to specified directory
#' @export
filter_empty_cells <- function(seurat_obj,
                               threshold = 100,
                               test.ambient = TRUE,
                               FDR = 0.001,
                               save_plots = TRUE,
                               plot_dir = "empty_cells_plots/",
                               plot_prefix = NULL,
                               verbose = TRUE,
                               BPPARAM = NULL,
                               ...) {
  
  # Load required packages
  required_packages <- c("DropletUtils", "ggplot2", "patchwork", "scales", "Seurat")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not available"))
    }
  }
  
  # Check if input is a list of Seurat objects or single object
  if (is.list(seurat_obj) && all(sapply(seurat_obj, function(x) inherits(x, "Seurat")))) {
    return(filter_empty_cells_list(seurat_obj, threshold, test.ambient, FDR, 
                                   save_plots, plot_dir, plot_prefix, verbose, 
                                   BPPARAM, ...))
  } else if (inherits(seurat_obj, "Seurat")) {
    return(filter_empty_cells_single(seurat_obj, threshold, test.ambient, FDR, 
                                     save_plots, plot_dir, plot_prefix, verbose, ...))
  } else {
    stop("Input must be a Seurat object or a list of Seurat objects")
  }
}

# =============================================================================
# SINGLE OBJECT PROCESSING FUNCTION
# =============================================================================

filter_empty_cells_single <- function(seurat_obj,
                                      threshold = 100,
                                      test.ambient = TRUE,
                                      FDR = 0.001,
                                      save_plots = TRUE,
                                      plot_dir = "empty_cells_plots/",
                                      plot_prefix = NULL,
                                      verbose = TRUE,
                                      ...) {
  
  if (verbose) {
    message("Starting empty droplet filtering...")
  }
  
  # Extract sample name from Seurat object
  sample_name <- seurat_obj@project.name
  
  # Get the original count matrix (this should be the RAW, unfiltered data)
  # Note: This assumes the Seurat object was created from raw data
  original_counts <- GetAssayData(seurat_obj, slot = "counts")
  
  if (verbose) {
    message(paste("Processing sample:", sample_name))
    message(paste("Original matrix:", nrow(original_counts), "genes x", ncol(original_counts), "cells"))
  }
  
  # Run emptyDrops analysis
  if (verbose) message("Running emptyDrops analysis...")
  
  empty_drops_results <- emptyDrops(
    original_counts, 
    lower = threshold, 
    test.ambient = test.ambient,
    ...
  )
  
  # Identify cells to keep
  keep_cells_idx <- which(empty_drops_results$FDR <= FDR)
  keep_cell_barcodes <- colnames(original_counts)[keep_cells_idx]
  
  # Calculate barcode ranks for all cells
  barcode_ranks_all <- barcodeRanks(original_counts)
  
  # Create metadata for ALL cells (before filtering)
  all_cell_metadata <- data.frame(
    cell_barcode = colnames(original_counts),
    total_UMI_raw = colSums(original_counts),
    genes_detected_raw = colSums(original_counts > 0),
    barcode_rank = barcode_ranks_all$rank,
    emptyDrops_FDR = empty_drops_results$FDR,
    emptyDrops_PValue = empty_drops_results$PValue,
    emptyDrops_LogProb = empty_drops_results$LogProb,
    cell_probability = 1 - empty_drops_results$FDR,  # Confidence it's a real cell
    passed_filter = 1:ncol(original_counts) %in% keep_cells_idx,
    stringsAsFactors = FALSE
  )
  
  # Handle NA values in cell_probability
  all_cell_metadata$cell_probability[is.na(all_cell_metadata$cell_probability)] <- 0
  
  # Create a NEW Seurat object starting with the RAW (unfiltered) data
  # This ensures both assays can have different numbers of cells
  options(Seurat.object.assay.version = "v5")
  filtered_seurat <- CreateSeuratObject(
    counts = original_counts,  # Start with ALL cells
    project = seurat_obj@project.name,
    min.cells = 0,  # Don't filter anything yet
    min.features = 0
  )
  
  # Add the RAW assay (same as the counts we just used)
  filtered_seurat[["RAW"]] <- CreateAssay5Object(counts = original_counts)
  
  # Add the RNA assay with ONLY the filtered cells
  filtered_counts <- original_counts[, keep_cells_idx, drop = FALSE]
  filtered_seurat[["RNA"]] <- CreateAssay5Object(counts = filtered_counts)
  
  # Set RNA as the default assay
  DefaultAssay(filtered_seurat) <- "RNA"
  
  # Add metadata for ALL cells (both kept and filtered out)
  # First, ensure rownames match cell barcodes
  rownames(all_cell_metadata) <- all_cell_metadata$cell_barcode
  
  # Add all metadata columns
  filtered_seurat$total_UMI_raw <- all_cell_metadata[colnames(filtered_seurat), "total_UMI_raw"]
  filtered_seurat$genes_detected_raw <- all_cell_metadata[colnames(filtered_seurat), "genes_detected_raw"]
  filtered_seurat$barcode_rank <- all_cell_metadata[colnames(filtered_seurat), "barcode_rank"]
  filtered_seurat$cell_probability <- all_cell_metadata[colnames(filtered_seurat), "cell_probability"]
  filtered_seurat$emptyDrops_FDR <- all_cell_metadata[colnames(filtered_seurat), "emptyDrops_FDR"]
  filtered_seurat$passed_emptyDrops <- all_cell_metadata[colnames(filtered_seurat), "passed_filter"]
  
  # Create summary statistics
  summary_stats <- list(
    sample_name = sample_name,
    cells_before = ncol(original_counts),
    cells_after = length(keep_cells_idx),
    cells_removed = ncol(original_counts) - length(keep_cells_idx),
    percent_kept = round(100 * length(keep_cells_idx) / ncol(original_counts), 2),
    threshold_used = threshold,
    FDR_threshold = FDR,
    genes_before = nrow(original_counts),
    genes_after = nrow(filtered_counts)
  )
  
  # Print summary
  if (verbose) {
    print_empty_drops_summary(summary_stats)
    message(paste("RAW assay contains", ncol(GetAssayData(filtered_seurat, assay = "RAW")), "cells (original data)"))
    message(paste("RNA assay contains", ncol(GetAssayData(filtered_seurat, assay = "RNA")), "cells (filtered data)"))
    message("Default assay set to 'RNA' (filtered cells)")
  }
  
  # Create plots
  plots_list <- create_empty_drops_plots(
    original_counts = original_counts,
    keep_cells_idx = keep_cells_idx,
    all_metadata = all_cell_metadata,
    barcode_ranks = barcode_ranks_all,
    sample_name = sample_name,
    threshold = threshold,
    save_plots = save_plots,
    plot_dir = plot_dir,
    plot_prefix = plot_prefix
  )
  
  if (verbose) {
    message("Empty cell filtering completed successfully!")
    message("Use GetAssayData(seurat_obj, assay='RAW') to access original counts")
    message("Use GetAssayData(seurat_obj, assay='RNA') to access filtered counts")
  }
  
  return(filtered_seurat)
}

# =============================================================================
# LIST PROCESSING FUNCTION
# =============================================================================

filter_empty_cells_list <- function(seurat_list,
                                    threshold = 100,
                                    test.ambient = TRUE,
                                    FDR = 0.001,
                                    save_plots = TRUE,
                                    plot_dir = "empty_cells_plots/",
                                    plot_prefix = NULL,
                                    verbose = TRUE,
                                    BPPARAM = NULL,
                                    ...) {
  
  if (verbose) {
    message(paste("Processing", length(seurat_list), "Seurat objects for empty cell filtering"))
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
    filtered_obj <- filter_empty_cells_single(
      seurat_obj = obj,
      threshold = threshold,
      test.ambient = test.ambient,
      FDR = FDR,
      save_plots = save_plots,
      plot_dir = plot_dir,
      plot_prefix = current_plot_prefix,
      verbose = verbose,
      ...
    )
    
    return(filtered_obj)
  }
  
  # Process objects (parallel or serial)
  if (!is.null(BPPARAM)) {
    # Parallel processing
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      warning("BiocParallel not available, falling back to serial processing")
      filtered_objects <- lapply(seq_along(seurat_list), process_one_object)
    } else {
      if (verbose) {
        message("Processing objects in parallel...")
      }
      
      filtered_objects <- tryCatch({
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
    filtered_objects <- lapply(seq_along(seurat_list), process_one_object)
  }
  
  # Restore names
  names(filtered_objects) <- sample_names
  
  if (verbose) {
    message(paste("Successfully filtered", length(filtered_objects), "Seurat objects"))
    
    # Print summary for all objects
    cat("\n=== Summary for All Samples ===\n")
    for (i in seq_along(filtered_objects)) {
      obj <- filtered_objects[[i]]
      name <- names(filtered_objects)[i]
      raw_cells <- ncol(GetAssayData(obj, assay = "RAW"))
      filtered_cells <- ncol(GetAssayData(obj, assay = "RNA"))
      percent_kept <- round(100 * filtered_cells / raw_cells, 2)
      
      cat(sprintf("%s: %d -> %d cells (%.1f%% kept)\n", 
                  name, raw_cells, filtered_cells, percent_kept))
    }
    cat("\n")
  }
  
  return(filtered_objects)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to print summary statistics
print_empty_drops_summary <- function(summary_stats) {
  cat("\n=== Empty Droplet Filtering Summary ===\n")
  cat("Sample:", summary_stats$sample_name, "\n")
  cat("Cells before filtering:", summary_stats$cells_before, "\n")
  cat("Cells after filtering:", summary_stats$cells_after, "\n")
  cat("Cells kept:", paste0(summary_stats$percent_kept, "%"), "\n")
  cat("Cells removed:", summary_stats$cells_removed, "\n")
  cat("Genes before filtering:", summary_stats$genes_before, "\n")
  cat("Genes after filtering:", summary_stats$genes_after, "\n")
  cat("UMI threshold used:", summary_stats$threshold_used, "\n")
  cat("FDR threshold used:", summary_stats$FDR_threshold, "\n\n")
}

# Enhanced function to create empty drops plots
create_empty_drops_plots <- function(original_counts, keep_cells_idx, all_metadata, 
                                     barcode_ranks, sample_name, threshold,
                                     save_plots = TRUE, plot_dir = "empty_cells_plots/",
                                     plot_prefix = NULL) {
  
  # Create plots directory if needed
  if (save_plots && !dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Plot 1: UMI distribution with filtering results
  p1 <- ggplot(all_metadata, aes(x = barcode_rank, y = total_UMI_raw, color = passed_filter)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue"), 
                       labels = c("FALSE" = "Filtered out", "TRUE" = "Kept"),
                       name = "Cell fate") +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "gray", size = 1) +
    labs(x = "Barcode rank", 
         y = "Total UMI count", 
         title = paste("Empty droplet filtering -", sample_name),
         subtitle = paste("Threshold:", threshold, "UMIs")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 2: Knee plot with filtering overlay
  knee_data <- data.frame(
    rank = barcode_ranks$rank,
    total = barcode_ranks$total,
    kept = seq_len(ncol(original_counts)) %in% keep_cells_idx
  )
  
  p2 <- ggplot(knee_data, aes(x = rank, y = total, color = kept)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue"), 
                       labels = c("FALSE" = "Filtered out", "TRUE" = "Kept"),
                       name = "Cell fate") +
    geom_hline(yintercept = metadata(barcode_ranks)$knee, color = "green", 
               linetype = "dashed", size = 1) +
    geom_hline(yintercept = metadata(barcode_ranks)$inflection, color = "orange", 
               linetype = "dashed", size = 1) +
    geom_hline(yintercept = threshold, color = "gray", linetype = "dashed", size = 1) +
    labs(x = "Barcode rank", 
         y = "Total UMI count",
         title = paste("Knee plot with filtering -", sample_name),
         subtitle = "Green: knee point, Orange: inflection point, Gray: threshold") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 3: Cell probability distribution
  p3 <- ggplot(all_metadata[all_metadata$passed_filter, ], aes(x = cell_probability)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = mean(all_metadata$cell_probability[all_metadata$passed_filter], na.rm = TRUE), 
               color = "red", linetype = "dashed", size = 1) +
    labs(x = "Cell Probability (1 - FDR)", 
         y = "Number of Cells",
         title = paste("Cell Probability Distribution -", sample_name),
         subtitle = paste("Mean probability:", 
                          round(mean(all_metadata$cell_probability[all_metadata$passed_filter], na.rm = TRUE), 3))) +
    theme_minimal()
  
  # Plot 4: Genes vs UMIs (for kept cells)
  kept_data <- all_metadata[all_metadata$passed_filter, ]
  p4 <- ggplot(kept_data, aes(x = total_UMI_raw, y = genes_detected_raw)) +
    geom_point(alpha = 0.6, size = 0.8, color = "blue") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = "Total UMI Count", 
         y = "Genes Detected",
         title = paste("Genes vs UMIs (Kept Cells) -", sample_name)) +
    theme_minimal()
  
  # Combine plots
  combined_plot <- (p1 | p2) / (p3 | p4) + 
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Save plot if requested
  if (save_plots) {
    if(is.null(plot_prefix)) {
      plot_filename <- file.path(plot_dir, paste0(sample_name, "_empty_cell_filtering.tiff"))
    } else {
      plot_filename <- file.path(plot_dir, paste0(plot_prefix, sample_name, "_empty_cell_filtering.tiff"))
    }
    
    ggsave(plot_filename, combined_plot, width = 14, height = 10, dpi = 300)
    message("Empty cell filtering plots saved to:", plot_filename)
  }
  
  # Return individual plots and combined plot
  plots_list <- list(
    umi_distribution = p1,
    knee_plot = p2,
    cell_probability = p3,
    genes_vs_umis = p4,
    combined = combined_plot,
    plot_filename = if(save_plots) plot_filename else NULL
  )
  
  return(plots_list)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Single Seurat object
# filtered_seurat <- filter_empty_cells(seurat_obj)

# List of Seurat objects (serial)
# seurat_list <- list(Sample1 = obj1, Sample2 = obj2, Sample3 = obj3)
# filtered_list <- filter_empty_cells(seurat_list)

# List of Seurat objects (parallel)
# library(BiocParallel)
# filtered_list <- filter_empty_cells(seurat_list, BPPARAM = MulticoreParam(workers = 4))

# Custom parameters for list
# filtered_list <- filter_empty_cells(
#   seurat_list,
#   threshold = 200,
#   FDR = 0.01,
#   plot_dir = "my_plots/",
#   plot_prefix = "experiment1_",
#   BPPARAM = MulticoreParam(workers = 2)
# )

# Access the different count matrices:
# raw_counts <- GetAssayData(filtered_seurat, assay = "RAW", slot = "counts")      # Original unfiltered
# filtered_counts <- GetAssayData(filtered_seurat, assay = "RNA", slot = "counts") # Filtered by emptyDrops

# For decontX background (empty droplets):
# raw_barcodes <- colnames(GetAssayData(filtered_seurat, assay = "RAW"))
# filtered_barcodes <- colnames(GetAssayData(filtered_seurat, assay = "RNA")) 
# empty_droplet_barcodes <- setdiff(raw_barcodes, filtered_barcodes)
# background_matrix <- GetAssayData(filtered_seurat, assay = "RAW")[, empty_droplet_barcodes]

# Check which cells passed filtering:
# table(filtered_seurat$passed_emptyDrops)  # TRUE = kept, FALSE = filtered out