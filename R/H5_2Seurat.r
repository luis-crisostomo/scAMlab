# Helper function to read CellBender H5 files and extract QC data
#' @export
read_cellbender_h5 <- function(filename) {
  if (!require(hdf5r)) stop("hdf5r package is required")
  
  h5_file <- H5File$new(filename, mode = "r")
  
  tryCatch({
    # CellBender structure: matrix/data, matrix/indices, matrix/indptr, matrix/shape
    # matrix/barcodes, matrix/features/name (or matrix/features/id)
    
    # Read the sparse matrix components
    data <- h5_file[["matrix/data"]][]
    indices <- h5_file[["matrix/indices"]][]
    indptr <- h5_file[["matrix/indptr"]][]
    shape <- h5_file[["matrix/shape"]][]
    
    # Read barcodes and features
    barcodes <- h5_file[["matrix/barcodes"]][]
    
    # Try to read feature names (gene symbols)
    features <- tryCatch({
      h5_file[["matrix/features/name"]][]
    }, error = function(e) {
      # Fallback to feature IDs if names don't exist
      tryCatch({
        h5_file[["matrix/features/id"]][]
      }, error = function(e2) {
        h5_file[["matrix/features"]][]
      })
    })
    
    # Convert to character if they're not already
    if (is.raw(barcodes)) barcodes <- rawToChar(barcodes)
    if (is.raw(features)) features <- rawToChar(features)
    
    # If they come as single strings, split them
    if (length(barcodes) == 1) {
      barcodes <- strsplit(barcodes, "\\0")[[1]]
      barcodes <- barcodes[barcodes != ""]
    }
    if (length(features) == 1) {
      features <- strsplit(features, "\\0")[[1]]
      features <- features[features != ""]
    }
    
    # Create sparse matrix
    mat <- Matrix::sparseMatrix(
      i = indices + 1,  # Convert from 0-based to 1-based indexing
      p = indptr,
      x = as.numeric(data),
      dims = shape,
      dimnames = list(features, barcodes)
    )
    
    # Extract QC data from CellBender
    qc_data <- list()
    
    # First, let's explore what's available in the file
    cat("Available groups in CellBender file:\n")
    tryCatch({
      all_names <- names(h5_file)
      print(all_names)
      
      # Check if droplet_latents exists
      if ("droplet_latents" %in% all_names) {
        droplet_names <- names(h5_file[["droplet_latents"]])
        cat("Available in droplet_latents:\n")
        print(droplet_names)
        
        # Read barcode indices for latents - this maps QC data to barcodes
        barcode_indices <- tryCatch({
          indices <- h5_file[["droplet_latents/barcode_indices_for_latents"]][]
          cat("Successfully read barcode_indices_for_latents with length:", length(indices), "\n")
          indices + 1  # Convert from 0-based to 1-based indexing
        }, error = function(e) {
          cat("Warning: Could not read barcode_indices_for_latents:", e$message, "\n")
          NULL
        })
        
        # Try to read cell probability
        qc_data$cell_probability <- tryCatch({
          cell_prob <- h5_file[["droplet_latents/cell_probability"]][]
          cat("Successfully read cell_probability with length:", length(cell_prob), "\n")
          list(values = cell_prob, indices = barcode_indices)
        }, error = function(e) {
          cat("Warning: Could not read cell_probability:", e$message, "\n")
          NULL
        })
        
        # Try to read other QC metrics if available
        if ("background_fraction" %in% droplet_names) {
          qc_data$background_fraction <- tryCatch({
            bg_frac <- h5_file[["droplet_latents/background_fraction"]][]
            list(values = bg_frac, indices = barcode_indices)
          }, error = function(e) { NULL })
        }
        
        if ("cell_size" %in% droplet_names) {
          qc_data$cell_size <- tryCatch({
            cell_sz <- h5_file[["droplet_latents/cell_size"]][]
            list(values = cell_sz, indices = barcode_indices)
          }, error = function(e) { NULL })
        }
        
        if ("droplet_efficiency" %in% droplet_names) {
          qc_data$droplet_efficiency <- tryCatch({
            drop_eff <- h5_file[["droplet_latents/droplet_efficiency"]][]
            list(values = drop_eff, indices = barcode_indices)
          }, error = function(e) { NULL })
        }
        
      } else {
        cat("droplet_latents group not found\n")
      }
      
      # Read global latents data
      if ("global_latents" %in% all_names) {
        global_names <- names(h5_file[["global_latents"]])
        cat("Available in global_latents:\n")
        print(global_names)
        
        # Try to read ambient expression
        if ("ambient_expression" %in% global_names) {
          qc_data$global_ambient_expression <- tryCatch({
            ambient_expr <- h5_file[["global_latents/ambient_expression"]][]
            cat("Successfully read global ambient_expression with length:", length(ambient_expr), "\n")
            ambient_expr
          }, error = function(e) {
            cat("Warning: Could not read global ambient_expression:", e$message, "\n")
            NULL
          })
        }
      } else {
        cat("global_latents group not found\n")
      }
    }, error = function(e) {
      cat("Error exploring file structure:", e$message, "\n")
    })
    
    h5_file$close_all()
    
    # Return both matrix and QC data
    return(list(matrix = mat, qc_data = qc_data, barcodes = barcodes))
    
  }, error = function(e) {
    h5_file$close_all()
    stop("Error reading CellBender H5 file structure: ", e$message)
  })
}

# Function to import CellBender output into Seurat object
#' @export
import_cellbender_to_seurat <- function(cellbender_h5_path, 
                                        raw_h5_path, 
                                        project_name = "CellBender_Analysis",
                                        min_cells = 3,
                                        min_features = 200) {
  
  # Load required libraries
  if (!require(Seurat)) stop("Seurat package is required")
  if (!require(hdf5r)) stop("hdf5r package is required")
  if (!require(Matrix)) stop("Matrix package is required")
  
  cat("Loading CellBender processed data...\n")
  
  # Import CellBender processed data (postCB.h5) as RNA assay
  # CellBender files have different structure, need to read manually
  tryCatch({
    # Try standard 10X format first
    cb_result <- tryCatch({
      cb_data <- Read10X_h5(filename = cellbender_h5_path, use.names = TRUE, unique.features = TRUE)
      list(matrix = cb_data, qc_data = list(), barcodes = colnames(cb_data))
    }, error = function(e) {
      # If that fails, read CellBender format manually
      cat("Standard 10X format failed, reading CellBender format...\n")
      read_cellbender_h5(cellbender_h5_path)
    })
    
    cb_data <- cb_result$matrix
    cb_qc <- cb_result$qc_data
    cb_barcodes <- cb_result$barcodes
    
    cat("✓ CellBender data loaded successfully\n")
    cat("  Dimensions:", nrow(cb_data), "genes x", ncol(cb_data), "cells\n")
    if (!is.null(cb_qc$cell_probability)) {
      cat("  QC data: cell_probability available for", length(cb_qc$cell_probability), "droplets\n")
    }
  }, error = function(e) {
    stop("Error loading CellBender data: ", e$message)
  })
  
  cat("Loading raw data...\n")
  
  # Import raw data as RAW assay
  tryCatch({
    raw_data <- Read10X_h5(filename = raw_h5_path, use.names = TRUE, unique.features = TRUE)
    cat("✓ Raw data loaded successfully\n")
    cat("  Dimensions:", nrow(raw_data), "genes x", ncol(raw_data), "cells\n")
  }, error = function(e) {
    stop("Error loading raw data: ", e$message)
  })
  
  cat("Creating Seurat object...\n")
  
  # Create Seurat object with CellBender data as the RNA assay
  seurat_obj <- CreateSeuratObject(
    counts = cb_data,
    project = project_name,
    assay = "RNA",
    min.cells = min_cells,
    min.features = min_features
  )
  
  cat("✓ Seurat object created with RNA assay\n")
  cat("  Final RNA assay dimensions:", 
      nrow(seurat_obj[["RNA"]]), "genes x", 
      ncol(seurat_obj[["RNA"]]), "cells\n")
  cat("  Built-in metadata columns:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
  
  # Get the barcodes that passed filtering in the RNA assay
  filtered_barcodes <- colnames(seurat_obj)
  
  # Subset raw data to match the filtered barcodes
  # This ensures both assays have the same cells
  common_barcodes <- intersect(filtered_barcodes, colnames(raw_data))
  
  if (length(common_barcodes) == 0) {
    stop("No common barcodes found between CellBender and raw data")
  }
  
  if (length(common_barcodes) < length(filtered_barcodes)) {
    warning("Some filtered barcodes not found in raw data. Using ", 
            length(common_barcodes), " common barcodes.")
    
    # Subset Seurat object to common barcodes
    seurat_obj <- subset(seurat_obj, cells = common_barcodes)
  }
  
  # Subset raw data to common barcodes and genes present in RNA assay
  raw_data_subset <- raw_data[rownames(seurat_obj), common_barcodes]
  
  cat("Adding RAW assay...\n")
  
  # Add raw data as RAW assay
  seurat_obj[["RAW"]] <- CreateAssayObject(counts = raw_data_subset)
  
  cat("✓ RAW assay added\n")
  cat("  RAW assay dimensions:", 
      nrow(seurat_obj[["RAW"]]), "genes x", 
      ncol(seurat_obj[["RAW"]]), "cells\n")
  
  # Add metadata about the analysis (preserve existing misc data)
  if (is.null(seurat_obj@misc$cellbender_info)) {
    seurat_obj@misc$cellbender_info <- list()
  }
  
  # Add/update cellbender info without overwriting
  seurat_obj@misc$cellbender_info$cellbender_file <- cellbender_h5_path
  seurat_obj@misc$cellbender_info$raw_file <- raw_h5_path
  seurat_obj@misc$cellbender_info$import_date <- Sys.time()
  seurat_obj@misc$cellbender_info$n_cells_raw <- ncol(raw_data)
  seurat_obj@misc$cellbender_info$n_cells_filtered <- ncol(seurat_obj)
  seurat_obj@misc$cellbender_info$n_genes_raw <- nrow(raw_data)
  seurat_obj@misc$cellbender_info$n_genes_filtered <- nrow(seurat_obj)
  
  # Add global ambient expression to misc if available
  if (!is.null(cb_qc$global_ambient_expression)) {
    seurat_obj@misc$cellbender_info$ambient_expression <- cb_qc$global_ambient_expression
    # If there were gene names, we could add them too for interpretation
    if (!is.null(cb_result$matrix)) {
      seurat_obj@misc$cellbender_info$ambient_expression_genes <- rownames(cb_result$matrix)
    }
    cat("✓ Added global ambient_expression to @misc slot\n")
    cat("  Ambient expression vector length:", length(cb_qc$global_ambient_expression), "\n")
  }
  
  # Add CellBender QC metadata to Seurat object
  # Match QC data to filtered barcodes using barcode indices
  
  # Helper function to map QC data to barcodes
  map_qc_to_barcodes <- function(qc_item, all_barcodes, final_barcodes) {
    if (is.null(qc_item) || is.null(qc_item$values) || is.null(qc_item$indices)) {
      return(NULL)
    }
    
    # Create a vector for all barcodes, initialized with NA
    full_qc <- rep(NA, length(all_barcodes))
    
    # Assign QC values to the correct barcode positions
    full_qc[qc_item$indices] <- qc_item$values
    
    # Name the vector with barcode names
    names(full_qc) <- all_barcodes
    
    # Return only the values for barcodes in the final Seurat object
    return(full_qc[final_barcodes])
  }
  
  if (!is.null(cb_qc$cell_probability)) {
    seurat_obj$cell_probability <- map_qc_to_barcodes(cb_qc$cell_probability, cb_barcodes, colnames(seurat_obj))
    if (!all(is.na(seurat_obj$cell_probability))) {
      cat("✓ Added cell_probability to metadata\n")
      cat("  Range:", round(range(seurat_obj$cell_probability, na.rm = TRUE), 4), "\n")
    }
  }
  
  if (!is.null(cb_qc$background_fraction)) {
    seurat_obj$background_fraction <- map_qc_to_barcodes(cb_qc$background_fraction, cb_barcodes, colnames(seurat_obj))
    if (!all(is.na(seurat_obj$background_fraction))) {
      cat("✓ Added background_fraction to metadata\n")
    }
  }
  
  if (!is.null(cb_qc$cell_size)) {
    seurat_obj$cell_size <- map_qc_to_barcodes(cb_qc$cell_size, cb_barcodes, colnames(seurat_obj))
    if (!all(is.na(seurat_obj$cell_size))) {
      cat("✓ Added cell_size to metadata\n")
    }
  }
  
  if (!is.null(cb_qc$droplet_efficiency)) {
    seurat_obj$droplet_efficiency <- map_qc_to_barcodes(cb_qc$droplet_efficiency, cb_barcodes, colnames(seurat_obj))
    if (!all(is.na(seurat_obj$droplet_efficiency))) {
      cat("✓ Added droplet_efficiency to metadata\n")
    }
  }
  
  # Calculate contamination metrics using Seurat's built-in nCount and nFeature
  # Seurat automatically creates nCount_RNA, nFeature_RNA, nCount_RAW, nFeature_RAW
  seurat_obj$contamination_ratio <- seurat_obj$nCount_RNA / seurat_obj$nCount_RAW
  seurat_obj$counts_removed <- seurat_obj$nCount_RAW - seurat_obj$nCount_RNA
  seurat_obj$features_removed <- seurat_obj$nFeature_RAW - seurat_obj$nFeature_RNA
  
  cat("\n=== Import Summary ===\n")
  cat("Project:", project_name, "\n")
  cat("Cells in final object:", ncol(seurat_obj), "\n")
  cat("Genes in final object:", nrow(seurat_obj), "\n")
  cat("Available assays:", names(seurat_obj@assays), "\n")
  cat("Default assay:", DefaultAssay(seurat_obj), "\n")
  
  # Print contamination summary
  cat("\n=== Contamination Summary ===\n")
  cat("Mean contamination ratio:", round(mean(seurat_obj$contamination_ratio), 4), "\n")
  cat("Median contamination ratio:", round(median(seurat_obj$contamination_ratio), 4), "\n")
  cat("Mean counts removed per cell:", round(mean(seurat_obj$counts_removed), 2), "\n")
  
  return(seurat_obj)
}

# Example usage function
example_usage <- function() {
  cat("Example usage:\n")
  cat("seurat_obj <- import_cellbender_to_seurat(\n")
  cat("  cellbender_h5_path = 'path/to/postCB.h5',\n")
  cat("  raw_h5_path = 'path/to/raw.h5',\n")
  cat("  project_name = 'MyProject',\n")
  cat("  min_cells = 3,\n")
  cat("  min_features = 200\n")
  cat(")\n\n")
  
  cat("# Access different assays:\n")
  cat("DefaultAssay(seurat_obj) <- 'RNA'  # Use CellBender-cleaned data\n")
  cat("DefaultAssay(seurat_obj) <- 'RAW'  # Use raw data\n\n")
  
  cat("# Compare contamination:\n")
  cat("VlnPlot(seurat_obj, features = c('n_counts_RNA', 'n_counts_RAW'))\n")
  cat("FeatureScatter(seurat_obj, feature1 = 'n_counts_RNA', feature2 = 'n_counts_RAW')\n")
}

# Helper function to compare assays
compare_assays <- function(seurat_obj) {
  if (!all(c("RNA", "RAW") %in% names(seurat_obj@assays))) {
    stop("Both RNA and RAW assays must be present")
  }
  
  cat("=== Assay Comparison ===\n")
  cat("RNA assay (CellBender cleaned):\n")
  cat("  Genes:", nrow(seurat_obj[["RNA"]]), "\n")
  cat("  Cells:", ncol(seurat_obj[["RNA"]]), "\n")
  cat("  Total counts:", sum(GetAssayData(seurat_obj, assay = "RNA", layer = "counts")), "\n")
  
  cat("RAW assay:\n")
  cat("  Genes:", nrow(seurat_obj[["RAW"]]), "\n")
  cat("  Cells:", ncol(seurat_obj[["RAW"]]), "\n")
  cat("  Total counts:", sum(GetAssayData(seurat_obj, assay = "RAW", layer = "counts")), "\n")
  
  cat("Difference:\n")
  rna_total <- sum(GetAssayData(seurat_obj, assay = "RNA", layer = "counts"))
  raw_total <- sum(GetAssayData(seurat_obj, assay = "RAW", layer = "counts"))
  
  cat("  Counts removed:", raw_total - rna_total, "\n")
  cat("  Percentage removed:", 
      round((1 - rna_total / raw_total) * 100, 2), "%\n")
}