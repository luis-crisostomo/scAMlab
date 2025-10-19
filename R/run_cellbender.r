library(Seurat)
library(reticulate)
library(Matrix)
library(hdf5r)

#' Helper function: Export Seurat object to h5 format for CellBender
#'
#' @param seurat_obj Seurat object containing raw count data
#' @param output_path Path for the output h5 file
#' @return Path to the created h5 file
#' @export
export_seurat_to_h5 <- function(seurat_obj, output_path) {
  
  cat("Exporting Seurat object to h5 format using DropletUtils...\n")
  
  # Get raw counts (compatible with both Seurat v4 and v5)
  if (!"RNA" %in% names(seurat_obj@assays)) {
    stop("RNA assay not found in Seurat object")
  }
  
  # Try to get raw counts using GetAssayData (works with both v4 and v5)
  tryCatch({
    raw_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
    
    # Check if we actually got count data
    if (is.null(raw_counts) || sum(raw_counts) == 0) {
      # Try the data slot as fallback
      raw_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
      if (!is.null(raw_counts) && sum(raw_counts) > 0) {
        warning("Using data slot instead of counts slot. Make sure these are raw counts!")
      } else {
        stop("No count data found in RNA assay")
      }
    }
    
    cat("Successfully extracted count matrix with", nrow(raw_counts), "features and", ncol(raw_counts), "cells\n")
    
  }, error = function(e) {
    # Fallback to direct slot access for older versions
    if (class(seurat_obj@assays$RNA)[1] == "Assay") {
      # Seurat v4 style
      if (!is.null(seurat_obj@assays$RNA@counts) && sum(seurat_obj@assays$RNA@counts) > 0) {
        raw_counts <- seurat_obj@assays$RNA@counts
      } else if (!is.null(seurat_obj@assays$RNA@data) && sum(seurat_obj@assays$RNA@data) > 0) {
        raw_counts <- seurat_obj@assays$RNA@data
        warning("Using data slot instead of counts slot. Make sure these are raw counts!")
      } else {
        stop("No count data found in RNA assay")
      }
    } else {
      stop("Could not access count data. Error: ", e$message)
    }
  })
  
  # Check for required package
  if (!requireNamespace("DropletUtils", quietly = TRUE)) {
    stop("DropletUtils package is required for h5 export. Please install it with: BiocManager::install('DropletUtils')")
  }
  
  # Remove existing file if it exists
  if (file.exists(output_path)) {
    file.remove(output_path)
    cat("Removed existing input file:", output_path, "\n")
  }
  
  tryCatch({
    cat("Writing h5 file with DropletUtils...\n")
    
    DropletUtils::write10xCounts(output_path, 
                                 raw_counts,
                                 barcodes = colnames(raw_counts),
                                 gene.id = rownames(raw_counts),
                                 gene.symbol = rownames(raw_counts),
                                 type = "HDF5",
                                 version = "3",
                                 overwrite = TRUE)
    
    cat("Successfully wrote h5 file to:", output_path, "\n")
    cat("File size:", round(file.size(output_path) / 1024^2, 2), "MB\n")
    
    # Verify the h5 file can be read back
    tryCatch({
      test_read <- Seurat::Read10X_h5(output_path)
      cat("H5 file validation successful:", nrow(test_read), "genes x", ncol(test_read), "cells\n")
    }, error = function(e) {
      warning("H5 file validation failed: ", e$message)
    })
    
    return(output_path)
    
  }, error = function(e) {
    stop("Failed to write h5 file with DropletUtils: ", e$message)
  })
}

#' Helper function: Run CellBender command
#'
#' @param input_h5_path Path to input h5 file
#' @param output_h5_path Path to output h5 file (without _filtered suffix)
#' @param output_dir Output directory
#' @param expected_cells Expected number of cells
#' @param total_droplets_included Total number of droplets to include
#' @param fpr False positive rate for cell detection
#' @param epochs Number of training epochs
#' @param cuda Use CUDA if available
#' @param learning_rate Learning rate
#' @param z_dim Latent dimension size
#' @param low_count_threshold Low count threshold
#' @param ... Additional CellBender parameters
#' @return List with exit_code and actual_output_path
#' @export
run_cellbender_command <- function(input_h5_path, output_h5_path, output_dir,
                                   expected_cells, total_droplets_included, fpr, epochs,
                                   cuda, learning_rate, z_dim, low_count_threshold, ...) {
  
  # Set environment variables for proper checkpoint handling
  Sys.setenv(CELLBENDER_CHECKPOINT_DIR = output_dir)
  cat("Set checkpoint directory to:", output_dir, "\n")
  
  # Prepare base CellBender command
  cellbender_cmd <- c(
    "remove-background",
    "--input", input_h5_path,
    "--output", output_h5_path,
    "--expected-cells", as.character(expected_cells),
    "--total-droplets-included", as.character(total_droplets_included),
    "--fpr", as.character(fpr),
    "--epochs", as.character(epochs),
    "--learning-rate", as.character(learning_rate),
    "--z-dim", as.character(z_dim),
    "--low-count-threshold", as.character(low_count_threshold)
  )
  
  # Add CUDA flag if requested and available
  if (cuda) {
    tryCatch({
      torch <- import("torch")
      if (torch$cuda$is_available()) {
        cellbender_cmd <- c(cellbender_cmd, "--cuda")
        cat("CUDA is available and will be used\n")
      } else {
        cat("CUDA requested but not available, running on CPU\n")
      }
    }, error = function(e) {
      cat("Could not check CUDA availability, running on CPU\n")
    })
  }
  
  # Add any additional parameters passed via ...
  extra_params <- list(...)
  if (length(extra_params) > 0) {
    for (param_name in names(extra_params)) {
      param_value <- extra_params[[param_name]]
      if (is.logical(param_value) && param_value) {
        # Boolean flags (like --debug)
        cellbender_cmd <- c(cellbender_cmd, paste0("--", gsub("_", "-", param_name)))
      } else if (!is.logical(param_value) || param_value != FALSE) {
        # Parameters with values
        cellbender_cmd <- c(cellbender_cmd, 
                            paste0("--", gsub("_", "-", param_name)), 
                            as.character(param_value))
      }
    }
  }
  
  cat("Running CellBender with command:\n")
  full_cmd <- paste("cellbender", paste(cellbender_cmd, collapse = " "))
  cat(full_cmd, "\n\n")
  
  # Execute the command
  result <- system(full_cmd, intern = FALSE)
  
  # Determine the actual output file (use the base name, not _filtered)
  expected_output <- output_h5_path  # Use the original path, not with _filtered suffix
  
  return(list(
    exit_code = result,
    expected_output_path = expected_output
  ))
}

#' Main function: Complete CellBender workflow
#'
#' @param seurat_obj Seurat object containing raw count data
#' @param output_dir Directory to save CellBender outputs
#' @param sample_name Name for the sample (will be used in output filenames)
#' @param expected_cells Expected number of cells (if NULL, will estimate from Seurat object)
#' @param total_droplets_included Total number of droplets to include (default: 25000)
#' @param fpr False positive rate for cell detection (default: 0.01)
#' @param epochs Number of training epochs (default: 150)
#' @param cuda Use CUDA if available (default: TRUE)
#' @param learning_rate Learning rate (default: 0.0001)
#' @param z_dim Latent dimension size (default: 64)
#' @param low_count_threshold Droplets with UMI counts below this will be completely excluded (default: 15)
#' @param min_cells Minimum number of cells expressing a feature for it to be kept (default: 5)
#' @param min_features Minimum number of features per cell (default: 200)
#' @param ... Additional CellBender parameters (e.g., debug = TRUE, checkpoint_mins = 30)
#' @return Seurat object with CellBender results (RNA = CellBender corrected, Raw = original)
#' @export
run_cellbender <- function(seurat_obj, 
                                      output_dir, 
                                      sample_name,
                                      expected_cells = NULL,
                                      total_droplets_included = 25000,
                                      fpr = 0.01,
                                      epochs = 150L,
                                      cuda = TRUE,
                                      learning_rate = 0.0001,
                                      z_dim = 64L,
                                      low_count_threshold = 15L,
                                      min_cells = 5,
                                      min_features = 200,
                                      ...) {
  
  cat("=== Starting Complete CellBender Analysis ===\n")
  cat("Sample:", sample_name, "\n")
  cat("Output directory:", output_dir, "\n\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Estimate expected cells if not provided
  if (is.null(expected_cells)) {
    expected_cells <- ncol(seurat_obj)
    cat("Estimated expected cells from Seurat object:", expected_cells, "\n")
  }
  
  # Prepare file paths
  input_h5_path <- file.path(output_dir, paste0(sample_name, "_raw.h5"))
  output_h5_path <- file.path(output_dir, paste0(sample_name, "_postCB.h5"))
  
  # Step 1: Export Seurat object to h5 format
  cat("\n=== Step 1: Exporting Seurat object to h5 format ===\n")
  input_path <- export_seurat_to_h5(seurat_obj, input_h5_path)
  
  # Step 2: Run CellBender
  cat("\n=== Step 2: Running CellBender ===\n")
  cb_result <- run_cellbender_command(
    input_h5_path = input_path,
    output_h5_path = output_h5_path,
    output_dir = output_dir,
    expected_cells = expected_cells,
    total_droplets_included = total_droplets_included,
    fpr = fpr,
    epochs = epochs,
    cuda = cuda,
    learning_rate = learning_rate,
    z_dim = z_dim,
    low_count_threshold = low_count_threshold,
    ...
  )
  
  if (cb_result$exit_code != 0) {
    stop("CellBender failed with exit code: ", cb_result$exit_code)
  }
  
  cat("CellBender completed successfully!\n")
  
  # Check if the expected output file exists
  if (!file.exists(cb_result$expected_output_path)) {
    # Look for other possible output files
    possible_outputs <- list.files(output_dir, pattern = paste0(sample_name, ".*\\.h5$"), full.names = TRUE)
    cat("Available h5 files in output directory:\n")
    print(possible_outputs)
    stop("Expected CellBender output file not found: ", cb_result$expected_output_path)
  }
  
  cat("CellBender output found at:", cb_result$expected_output_path, "\n")
  
  # Step 3: Create integrated Seurat object using your working function
  cat("\n=== Step 3: Creating integrated Seurat object ===\n")
  
  # This section requires function "import_cellbender_to_seurat"
  
  integrated_seurat <- import_cellbender_to_seurat(
    cellbender_h5_path = cb_result$expected_output_path, 
    raw_h5_path = input_path, 
    project_name = sample_name,
    min_cells = min_cells,
    min_features = min_features
  )
  
  # Add metadata about the CellBender run
  integrated_seurat@misc$cellbender_params <- list(
    sample_name = sample_name,
    expected_cells = expected_cells,
    total_droplets_included = total_droplets_included,
    fpr = fpr,
    epochs = epochs,
    learning_rate = learning_rate,
    z_dim = z_dim,
    low_count_threshold = low_count_threshold,
    min_cells = min_cells,
    min_features = min_features,
    output_dir = output_dir,
    input_file = input_path,
    output_file = cb_result$expected_output_path,
    run_date = Sys.time()
  )
  
  # Add comparison summary
  cat("\n=== Comparison Summary ===\n")
  cat("Original Seurat object: ", ncol(seurat_obj), "cells,", nrow(seurat_obj), "features\n")
  cat("CellBender filtered (RNA assay):", ncol(integrated_seurat@assays$RNA), "cells,", nrow(integrated_seurat@assays$RNA), "features\n")
  cat("Raw counts (RAW assay):", ncol(integrated_seurat@assays$RAW), "cells,", nrow(integrated_seurat@assays$RAW), "features\n")
  
  cat("\n=== Complete CellBender Analysis Complete ===\n")
  cat("Returning Seurat object with:\n")
  cat("- RNA assay: CellBender corrected counts\n")
  cat("- RAW assay: Original raw counts\n")
  
  return(integrated_seurat)
  invisible(file.remove("ckpt.tar.gz"))
}

#' Convenience wrapper function for complete workflow with comparison
#'
#' @param seurat_obj Input Seurat object
#' @param output_dir Output directory
#' @param sample_name Sample name
#' @param ... Additional parameters passed to run_cellbender
#' @return List containing original Seurat object and CellBender-processed object
#' @export
run_cellbender_with_comparison <- function(seurat_obj, output_dir, sample_name, ...) {
  
  cat("Starting complete CellBender workflow for sample:", sample_name, "\n")
  
  # Run integrated CellBender workflow
  cellbender_seurat <- run_cellbender(
    seurat_obj = seurat_obj,
    output_dir = output_dir,
    sample_name = sample_name,
    ...
  )
  
  return(list(
    original = seurat_obj,
    integrated = cellbender_seurat,
    output_dir = output_dir
  ))
  invisible(file.remove("ckpt.tar.gz"))
}
  