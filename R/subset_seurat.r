#' Run scDblFinder for Doublet Detection in Single-Cell RNA-seq Data
#'
#' @description
#' Detects doublets (artificial cells formed by two or more cells captured
#' together) in single-cell RNA-seq data using scDblFinder. The function
#' converts Seurat objects to SingleCellExperiment format, optionally performs
#' fast clustering for improved doublet detection, and adds doublet scores and
#' classifications back to the Seurat object. Supports both single Seurat
#' objects and lists of Seurat objects with flexible parameter specification.
#'
#' @param seurat_obj A Seurat object or a list of Seurat objects to process
#'   for doublet detection
#' @param clusters Logical. If TRUE, runs fast clustering before doublet
#'   detection to improve accuracy by identifying cell populations. If FALSE,
#'   performs doublet detection without clustering. Can be a single value
#'   (applied to all objects) or a vector matching the length of seurat_obj
#'   when processing a list. Default: TRUE
#' @param assay Character string. Name of the assay to use for doublet detection.
#'   Can be a single value (applied to all objects) or a vector matching the
#'   length of seurat_obj when processing a list. Default: "RNA"
#' @param ... Additional arguments passed to \code{\link[scDblFinder]{scDblFinder}}.
#'   Common arguments include:
#'   \itemize{
#'     \item \strong{dbr}: Expected doublet rate (NULL for automatic estimation)
#'     \item \strong{dbr.sd}: Standard deviation of doublet rate (default: 0.015)
#'     \item \strong{nfeatures}: Number of top features to use (default: 1000)
#'     \item \strong{dims}: Number of PCA dimensions to use (default: 20)
#'     \item \strong{k}: Number of nearest neighbors (default: 10)
#'     \item \strong{BPPARAM}: BiocParallelParam for parallel processing
#'   }
#'   When processing a list, can be single values (applied to all) or vectors/lists
#'   matching the number of objects.
#'
#' @return
#' \itemize{
#'   \item For single Seurat object: Returns the input Seurat object with added
#'         metadata columns:
#'     \itemize{
#'       \item \strong{doublet_scores}: Numeric doublet scores (0-1, higher = more likely doublet)
#'       \item \strong{doublet_class}: Factor classification ("singlet" or "doublet")
#'       \item \strong{fastcluster_clusters}: Factor cluster assignments (if clusters = TRUE)
#'     }
#'   \item For list of Seurat objects: Returns a named list of processed Seurat
#'         objects with the same structure as above
#' }
#'
#' @details
#' **Doublet Detection Overview:**
#'
#' Doublets are a common technical artifact in droplet-based scRNA-seq where
#' two or more cells are captured in the same droplet, creating artificial
#' transcriptional profiles. scDblFinder uses:
#' \enumerate{
#'   \item Simulation of artificial doublets from the data
#'   \item Machine learning to distinguish real doublets from singlets
#'   \item Cluster-aware detection to account for cell type composition
#' }
#'
#' **Why Use Clustering?**
#'
#' Setting \code{clusters = TRUE} (recommended) provides several advantages:
#' \itemize{
#'   \item \strong{Improved accuracy}: Doublet detection considers cell type composition
#'   \item \strong{Reduced false positives}: Prevents rare cell types from being
#'         flagged as doublets
#'   \item \strong{Better doublet simulation}: Creates more realistic artificial
#'         doublets based on actual cell populations
#'   \item \strong{Cluster information}: Adds fastcluster_clusters to metadata
#'         for QC visualization
#' }
#'
#' Set \code{clusters = FALSE} only if:
#' \itemize{
#'   \item Dataset is very small (< 500 cells)
#'   \item All cells are expected to be from one homogeneous population
#'   \item You want faster processing and already have cluster information
#' }
#'
#' **Expected Doublet Rates:**
#'
#' Typical doublet rates by technology (if dbr not specified, scDblFinder estimates):
#' \itemize{
#'   \item 10X Chromium: ~1% per 1,000 cells loaded (e.g., 8% at 8,000 cells)
#'   \item Drop-seq: ~2% per 1,000 cells
#'   \item inDrops: ~2% per 1,000 cells
#' }
#'
#' The actual rate depends on cell loading concentration and capture efficiency.
#'
#' **Processing Multiple Objects:**
#'
#' When processing a list of Seurat objects, parameters can be specified as:
#' \itemize{
#'   \item \strong{Single values}: Applied to all objects (e.g., clusters = TRUE)
#'   \item \strong{Vectors}: Different value for each object (e.g., assay = c("RNA", "SCT", "RNA"))
#'   \item \strong{Lists}: For complex parameters like BPPARAM
#' }
#'
#' All parameters are validated to ensure correct length matching the number of objects.
#'
#' **Memory Management:**
#'
#' The function includes automatic memory management:
#' \itemize{
#'   \item Garbage collection after conversion to SingleCellExperiment
#'   \item Removal of temporary SCE object after processing
#'   \item Garbage collection between objects when processing lists
#' }
#'
#' **Interpreting Results:**
#'
#' \strong{doublet_scores}:
#' \itemize{
#'   \item Range: 0 to 1
#'   \item Higher scores indicate higher probability of being a doublet
#'   \item Typical threshold: 0.5 (but examine distribution in your data)
#' }
#'
#' \strong{doublet_class}:
#' \itemize{
#'   \item "singlet": Predicted single cell
#'   \item "doublet": Predicted doublet (typically scores > 0.5)
#'   \item Based on automatic threshold selection by scDblFinder
#' }
#'
#' @section Workflow Integration:
#' Doublet detection should typically be performed:
#' \enumerate{
#'   \item After basic QC (gene/cell filtering)
#'   \item Before normalization and integration
#'   \item Before clustering for downstream analysis
#' }
#'
#' Recommended workflow:
#' \preformatted{
#' # 1. Basic filtering
#' seurat_obj <- lowpass_gene_filter(seurat_obj, min_cells = 3)
#'
#' # 2. Doublet detection
#' seurat_obj <- run_scDblFinder(seurat_obj, clusters = TRUE)
#'
#' # 3. Visualize doublets
#' plot_QC_UMAP(seurat_obj)  # Check doublet distribution
#'
#' # 4. Remove doublets
#' seurat_obj <- subset(seurat_obj, subset = doublet_class == "singlet")
#'
#' # 5. Continue with normalization and analysis
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' seurat_obj <- run_scDblFinder(seurat_obj)
#'
#' # Check results
#' table(seurat_obj$doublet_class)
#' summary(seurat_obj$doublet_scores)
#'
#' # Visualize doublet scores
#' FeaturePlot(seurat_obj, features = "doublet_scores")
#'
#' # Without clustering (faster but less accurate)
#' seurat_obj <- run_scDblFinder(seurat_obj, clusters = FALSE)
#'
#' # Specify expected doublet rate (e.g., 5%)
#' seurat_obj <- run_scDblFinder(
#'   seurat_obj = seurat_obj,
#'   dbr = 0.05
#' )
#'
#' # Use different assay
#' seurat_obj <- run_scDblFinder(
#'   seurat_obj = seurat_obj,
#'   assay = "SCT"
#' )
#'
#' # Custom parameters for more sensitive detection
#' seurat_obj <- run_scDblFinder(
#'   seurat_obj = seurat_obj,
#'   nfeatures = 2000,  # Use more features
#'   dims = 30,         # Use more PCA dimensions
#'   k = 20             # Use more neighbors
#' )
#'
#' # Parallel processing with BiocParallel
#' library(BiocParallel)
#' seurat_obj <- run_scDblFinder(
#'   seurat_obj = seurat_obj,
#'   BPPARAM = MulticoreParam(workers = 4)
#' )
#'
#' # Process multiple samples
#' seurat_list <- list(
#'   Sample1 = seurat_obj1,
#'   Sample2 = seurat_obj2,
#'   Sample3 = seurat_obj3
#' )
#'
#' # Same parameters for all
#' results <- run_scDblFinder(seurat_list, clusters = TRUE)
#'
#' # Different parameters for each sample
#' results <- run_scDblFinder(
#'   seurat_obj = seurat_list,
#'   clusters = c(TRUE, TRUE, FALSE),  # No clustering for sample 3
#'   assay = c("RNA", "RNA", "SCT"),   # Different assay for sample 3
#'   dbr = c(0.05, 0.08, 0.06)         # Different doublet rates
#' )
#'
#' # Remove doublets from all samples
#' singlets_list <- lapply(results, function(obj) {
#'   subset(obj, subset = doublet_class == "singlet")
#' })
#'
#' # Check doublet rates across samples
#' doublet_summary <- sapply(results, function(obj) {
#'   sum(obj$doublet_class == "doublet") / ncol(obj)
#' })
#' print(doublet_summary)
#' }
#'
#' @references
#' Germain, P. L., Lun, A., Garcia Meixide, C., Macnair, W., & Robinson, M. D. (2021).
#' Doublet identification in single-cell sequencing data using scDblFinder.
#' F1000Research, 10, 979. DOI: https://doi.org/10.12688/f1000research.73600.2
#'
#' scDblFinder documentation: https://bioconductor.org/packages/scDblFinder/
#'
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom scDblFinder scDblFinder fastcluster
#' @importFrom SingleCellExperiment colData
#'
#' @seealso
#' \code{\link[scDblFinder]{scDblFinder}} for the underlying doublet detection algorithm
#' \code{\link[scDblFinder]{fastcluster}} for fast clustering implementation
#' \code{\link{plot_QC_UMAP}} for visualizing doublet detection results
#' \code{\link{run_babraham_qc}} for comprehensive QC before doublet detection
#'
#' @export
subset_seurat <- function(seurat_obj, subset_criteria, verbose = TRUE, sample = NULL,
                          assay = NULL, param = NULL) {

  # Function to process a single Seurat object
  process_single_seurat <- function(obj, criteria, obj_name = NULL, target_assay = NULL,
                                    be_verbose = TRUE) {

    # Validate input
    if (!inherits(obj, "Seurat")) {
      stop(paste("Object", ifelse(is.null(obj_name), "", paste0("'", obj_name, "'")),
                 "is not a Seurat object"))
    }

    if (!is.list(criteria) || length(criteria) == 0) {
      stop("subset_criteria must be a non-empty list")
    }

    # Check Seurat version and encourage v5 usage
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      if (be_verbose) warning("SeuratObject package not available. Cannot check Seurat version.")
    } else {
      # Check if object has v5 assay structure
      has_v5_assays <- any(sapply(obj@assays, function(x) inherits(x, "Assay5")))
      if (!has_v5_assays && be_verbose) {
        message("This function is optimized for Seurat v5 objects with Assay5 structure.")
        message("Consider updating your object using: UpdateSeuratObject() or JoinLayers()")
        message("Proceeding with standard subsetting...")
      }
    }

    # Handle assay selection
    if (is.null(target_assay)) {
      current_assay <- DefaultAssay(obj)
    } else {
      if (!target_assay %in% names(obj@assays)) {
        stop(paste("Assay '", target_assay, "' not found in",
                   ifelse(is.null(obj_name), "Seurat object", paste0("object '", obj_name, "'")),
                   "\nAvailable assays:", paste(names(obj@assays), collapse = ", ")))
      }
      current_assay <- target_assay
    }

    # Get available metadata columns
    available_cols <- colnames(obj@meta.data)

    # Store original cell count
    original_cells <- ncol(obj)

    # Get sample name
    if (is.null(obj_name)) {
      obj_name <- deparse(substitute(obj))
    }

    if (be_verbose) {
      cat("Starting QC filtering...\n")
      cat("Sample:", obj_name,"\n")
      cat("Target assay:", current_assay, "\n")
      cat("Original number of cells:", original_cells, "\n")
      cat("Available metadata columns:", paste(available_cols, collapse = ", "), "\n\n")
    }

    # Build logical vectors for each filtering condition
    filter_conditions <- rep(TRUE, ncol(obj))
    applied_filters <- c()

    for (param_name in names(criteria)) {

      # Check if parameter exists in metadata
      if (!param_name %in% available_cols) {
        warning(paste("Parameter '", param_name, "' not found in", obj_name, "metadata.",
                      "\nAvailable columns:", paste(available_cols, collapse = ", ")))
        next
      }

      param_conditions <- criteria[[param_name]]
      param_data <- obj@meta.data[[param_name]]

      # Handle different condition types
      if ("min" %in% names(param_conditions)) {
        condition <- param_data >= param_conditions$min
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters, paste(param_name, ">=", param_conditions$min))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, ">=", param_conditions$min,
              "- Cells passing:", cells_passing, "\n")
        }
      }

      if ("max" %in% names(param_conditions)) {
        condition <- param_data <= param_conditions$max
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters, paste(param_name, "<=", param_conditions$max))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, "<=", param_conditions$max,
              "- Cells passing:", cells_passing, "\n")
        }
      }

      if ("greater_than" %in% names(param_conditions)) {
        condition <- param_data > param_conditions$greater_than
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters, paste(param_name, ">", param_conditions$greater_than))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, ">", param_conditions$greater_than,
              "- Cells passing:", cells_passing, "\n")
        }
      }

      if ("less_than" %in% names(param_conditions)) {
        condition <- param_data < param_conditions$less_than
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters, paste(param_name, "<", param_conditions$less_than))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, "<", param_conditions$less_than,
              "- Cells passing:", cells_passing, "\n")
        }
      }

      if ("equals" %in% names(param_conditions)) {
        condition <- param_data == param_conditions$equals
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters, paste(param_name, "==", param_conditions$equals))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, "==", param_conditions$equals,
              "- Cells passing:", cells_passing, "\n")
        }
      }

      if ("not_equals" %in% names(param_conditions)) {
        condition <- param_data != param_conditions$not_equals
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters, paste(param_name, "!=", param_conditions$not_equals))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, "!=", param_conditions$not_equals,
              "- Cells passing:", cells_passing, "\n")
        }
      }

      # Handle range (between two values)
      if ("range" %in% names(param_conditions)) {
        if (length(param_conditions$range) != 2) {
          warning(paste("Range for", param_name, "must have exactly 2 values. Skipping."))
          next
        }
        condition <- param_data >= param_conditions$range[1] &
          param_data <= param_conditions$range[2]
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters,
                             paste(param_name, "between", param_conditions$range[1],
                                   "and", param_conditions$range[2]))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, "between", param_conditions$range[1],
              "and", param_conditions$range[2], "- Cells passing:", cells_passing, "\n")
        }
      }

      # Handle factor filtering - keep only specified levels
      if ("keep_levels" %in% names(param_conditions)) {
        levels_to_keep <- param_conditions$keep_levels
        if (!is.vector(levels_to_keep)) {
          warning(paste("keep_levels for", param_name, "must be a vector. Skipping."))
          next
        }
        condition <- param_data %in% levels_to_keep
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters,
                             paste(param_name, "in", paste(levels_to_keep, collapse = ", ")))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, "keep levels:", paste(levels_to_keep, collapse = ", "),
              "- Cells passing:", cells_passing, "\n")
        }
      }

      # Handle factor filtering - exclude specified levels
      if ("exclude_levels" %in% names(param_conditions)) {
        levels_to_exclude <- param_conditions$exclude_levels
        if (!is.vector(levels_to_exclude)) {
          warning(paste("exclude_levels for", param_name, "must be a vector. Skipping."))
          next
        }
        condition <- !param_data %in% levels_to_exclude
        filter_conditions <- filter_conditions & condition
        applied_filters <- c(applied_filters,
                             paste(param_name, "exclude", paste(levels_to_exclude, collapse = ", ")))
        if (be_verbose) {
          cells_passing <- sum(condition)
          cat("Filter:", param_name, "exclude levels:", paste(levels_to_exclude, collapse = ", "),
              "- Cells passing:", cells_passing, "\n")
        }
      }
    }

    # Check if any valid filters were applied
    if (length(applied_filters) == 0) {
      warning(paste("No valid filtering parameters found for", obj_name, ". Returning original object."))
      return(obj)
    }

    if (be_verbose) {
      cat("\nApplied filters:\n")
      for (filter_desc in applied_filters) {
        cat("-", filter_desc, "\n")
      }
    }

    # Apply filtering by subsetting cells
    tryCatch({
      # Get cell names that pass all filters
      cells_to_keep <- colnames(obj)[filter_conditions]

      # Subset the Seurat object with proper assay handling
      if (!is.null(target_assay)) {
        # Set the target assay as active, then subset
        DefaultAssay(obj) <- current_assay
        filtered_obj <- subset(obj, cells = cells_to_keep)
      } else {
        # Standard subsetting with currently active assay
        filtered_obj <- subset(obj, cells = cells_to_keep)
      }

      final_cells <- ncol(filtered_obj)
      removed_cells <- original_cells - final_cells
      percent_remaining <- round((final_cells / original_cells) * 100, 2)

      if (be_verbose) {
        cat("\nFiltering completed successfully!\n")
        cat("Cells remaining:", final_cells, "\n")
        cat("Cells removed:", removed_cells, "\n")
        cat("Percentage remaining:", percent_remaining, "%\n")
        cat("Active assay in result:", DefaultAssay(filtered_obj), "\n")
        cat("All layers within the assay have been subset consistently.\n")
        cat("========\n")
      }

      return(filtered_obj)

    }, error = function(e) {
      stop(paste("Error during filtering of", obj_name, ":", e$message,
                 "\nPlease check that parameter names exist in your Seurat object metadata.",
                 "\nUse list_metadata() to see available columns.",
                 "\nFor Seurat v5 objects, all layers within the assay will be subset together."))
    })
  }

  # Helper function to check if a list represents filtering criteria
  is_criteria_list <- function(x) {
    # A criteria list should have named elements where each element contains
    # filtering conditions (min, max, equals, etc.)
    if (!is.list(x) || length(x) == 0) return(FALSE)
    if (is.null(names(x))) return(FALSE)

    # Check if any element contains typical filtering conditions
    condition_names <- c("min", "max", "greater_than", "less_than", "equals",
                         "not_equals", "range", "keep_levels", "exclude_levels")

    # At least one element should contain filtering conditions
    has_conditions <- any(sapply(x, function(element) {
      if (!is.list(element)) return(FALSE)
      any(names(element) %in% condition_names)
    }))

    return(has_conditions)
  }

  # Main function logic starts here

  # Check if input is a list of Seurat objects
  if (is.list(seurat_obj) && !inherits(seurat_obj, "Seurat")) {
    # Handle list of Seurat objects
    if (length(seurat_obj) == 0) {
      stop("Input list is empty")
    }

    # Get names for the list elements
    obj_names <- names(seurat_obj)
    if (is.null(obj_names)) {
      obj_names <- paste0("Object_", seq_along(seurat_obj))
      names(seurat_obj) <- obj_names
    }

    # Handle criteria - can be single list or list of lists
    if (is.list(subset_criteria) && length(subset_criteria) > 0) {
      # Check if it's a single criteria list (to be applied to all objects)
      # or a list of criteria lists (different criteria for each object)

      if (is_criteria_list(subset_criteria)) {
        # Single criteria list for all objects
        criteria_list <- rep(list(subset_criteria), length(seurat_obj))
        if (verbose) {
          cat("Applying the same filtering criteria to all", length(seurat_obj), "objects.\n")
        }
      } else if (all(sapply(subset_criteria, is_criteria_list))) {
        # List of criteria lists (different criteria for each object)
        if (length(subset_criteria) != length(seurat_obj)) {
          stop("When providing a list of criteria lists, it must have the same length as the list of Seurat objects")
        }
        criteria_list <- subset_criteria
        if (verbose) {
          cat("Applying different filtering criteria to each of", length(seurat_obj), "objects.\n")
        }
      } else {
        stop("subset_criteria must be either:\n",
             "1) A single criteria list (applied to all objects), or\n",
             "2) A list of criteria lists (one for each object)\n",
             "Each criteria list should contain filtering conditions like min, max, equals, etc.")
      }
    } else {
      stop("subset_criteria must be a list or list of lists")
    }

    # Process objects
    if (!is.null(param)) {
      # Parallel processing using BiocParallel param
      if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        warning("BiocParallel package not available. Proceeding with sequential processing.")
      } else {
        # Get number of workers from BiocParallel param object
        n_workers <- tryCatch({
          BiocParallel::bpnworkers(param)
        }, error = function(e) {
          warning("Could not extract number of workers from param object. Proceeding with sequential processing.")
          return(NULL)
        })

        if (!is.null(n_workers) && n_workers > 1) {
          if (verbose) {
            cat("Using BiocParallel with", n_workers, "workers...\n\n")
          }

          # Ensure Seurat is loaded on all workers
          tryCatch({
            BiocParallel::bplapply(1:n_workers, function(x) {
              library(Seurat)
              if (requireNamespace("SeuratObject", quietly = TRUE)) {
                library(SeuratObject)
              }
              return(TRUE)
            }, BPPARAM = param)
          }, error = function(e) {
            warning("Failed to load Seurat on workers: ", e$message)
          })

          # Process in parallel using BiocParallel
          results <- tryCatch({
            BiocParallel::bplapply(seq_along(seurat_obj), function(i) {
              # Load libraries again to be safe
              if (!requireNamespace("Seurat", quietly = TRUE)) {
                stop("Seurat package not available on worker")
              }
              library(Seurat)
              if (requireNamespace("SeuratObject", quietly = TRUE)) {
                library(SeuratObject)
              }

              process_single_seurat(
                obj = seurat_obj[[i]],
                criteria = criteria_list[[i]],
                obj_name = obj_names[i],
                target_assay = assay,
                be_verbose = FALSE  # Turn off verbose for parallel
              )
            }, BPPARAM = param)
          }, error = function(e) {
            warning("Parallel processing failed: ", e$message, ". Falling back to sequential processing.")
            return(NULL)
          })

          if (!is.null(results)) {
            names(results) <- obj_names

            if (verbose) {
              cat("\n==========================================\n")
              cat("All objects processed successfully!\n")
              cat("Processed", length(results), "objects.\n")
            }

            return(results)
          }
        }
      }
    }

    # Sequential processing (fallback or when param not provided)
    results <- list()
    for (i in seq_along(seurat_obj)) {
      if (verbose) {
        cat("\n=== Processing:", obj_names[i], "===\n")
      }
      results[[obj_names[i]]] <- process_single_seurat(
        obj = seurat_obj[[i]],
        criteria = criteria_list[[i]],
        obj_name = obj_names[i],
        target_assay = assay,
        be_verbose = verbose
      )
    }

    names(results) <- obj_names

    if (verbose) {
      cat("\n==========================================\n")
      cat("All objects processed successfully!\n")
      cat("Processed", length(results), "objects.\n")
    }

    return(results)

  } else {
    # Handle single Seurat object
    return(process_single_seurat(
      obj = seurat_obj,
      criteria = subset_criteria,
      obj_name = sample,
      target_assay = assay,
      be_verbose = verbose
    ))
  }
}
