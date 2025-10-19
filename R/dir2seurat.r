# =============================================================================
# CLEAN dir2seurat FUNCTION - DATA IMPORT ONLY
# =============================================================================

#' Import single-cell RNA-seq data from directories into Seurat objects
#' 
#' This function loads raw count matrices from various single-cell formats
#' (10X, Parse Biosciences, HDF5) and creates Seurat objects. No filtering
#' or preprocessing is applied - this is pure data import.
#' 
#' @param data_dir Character vector of directory paths containing sc-RNA data
#' @param data_type Character, data format type ("auto", "tenx", "parse", "h5"). 
#'   Default "auto" attempts automatic detection.
#' @param sample_name Character vector of sample names. If NULL, inferred from directory names.
#' @param min_cells Minimum number of cells expressing a gene (passed to CreateSeuratObject)
#' @param min_features Minimum number of features per cell (passed to CreateSeuratObject)
#' @param project Project name for Seurat object
#' @param BPPARAM BiocParallel parameter for parallel processing (optional)
#' @param ... Additional arguments passed to data loading functions
#' 
#' @return Single Seurat object (if one directory) or named list of Seurat objects (multiple directories)
#' @export
dir2Seurat <- function(data_dir, 
                       data_type = "auto", 
                       sample_name = NULL,
                       min_cells = 3,
                       min_features = 200, 
                       project = "SingleCellProject",
                       BPPARAM = NULL,
                       ...) {
  
  # Load required packages upfront
  required_packages <- c("Seurat")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not available"))
    }
  }
  
  # Validate inputs
  if (!all(dir.exists(data_dir))) {
    missing_dirs <- data_dir[!dir.exists(data_dir)]
    stop(paste("Directory(ies) not found:", paste(missing_dirs, collapse = ", ")))
  }
  
  if (!is.null(sample_name) && length(sample_name) != length(data_dir)) {
    stop("Length of sample_name must match length of data_dir")
  }
  
  message(paste("Loading data from", length(data_dir), "directory(ies)"))
  
  # Internal function to process one directory
  process_one_directory <- function(dir_path, sample_name_single = NULL, counts_data_preloaded = NULL) {
  
    # Determine sample name
    if (is.null(sample_name_single)) {
      sample_name_single <- basename(normalizePath(dir_path))
    }
    
    message(paste("Processing sample:", sample_name_single))
    
    # Use preloaded counts or load them
    if (!is.null(counts_data_preloaded)) {
      counts_data <- counts_data_preloaded
    } else {
      # Load counts matrix
      counts_data <- get_count_matrix(dir_path, data_type = data_type, 
                                      sample_name = sample_name_single, ...)
    }
    
    # Create basic Seurat object
    options(Seurat.object.assay.version = "v5")
    
    seurat_obj <- CreateSeuratObject(
      counts = counts_data, 
      project = ifelse(length(data_dir) == 1, project, sample_name_single),
      min.cells = min_cells,
      min.features = min_features
    )
    
    message(paste("Created Seurat object for", sample_name_single, 
                  "with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes"))
    
    return(seurat_obj)
  }
  
  # Process single directory
  if (length(data_dir) == 1) {
    result <- process_one_directory(data_dir, sample_name)
    return(result)
  }
  
  # Process multiple directories
  if (!is.null(BPPARAM)) {
    # Parallel processing
    
    # Pre-load all count matrices serially first (file I/O can be tricky in parallel)
    message("Pre-loading count matrices...")
    all_counts <- lapply(seq_along(data_dir), function(i) {
      sample_name_single <- if (is.null(sample_name)) NULL else sample_name[i]
      get_count_matrix(data_dir[i], data_type = data_type, sample_name = sample_name_single, ...)
    })
    
    # Then process in parallel with preloaded data
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      warning("BiocParallel not available, falling back to serial processing")
      objs <- lapply(seq_along(data_dir), function(i) {
        process_one_directory(data_dir[i], 
                              if (is.null(sample_name)) NULL else sample_name[i],
                              all_counts[[i]])
      })
    } else {
      objs <- tryCatch({
        BiocParallel::bplapply(seq_along(data_dir), function(i) {
          process_one_directory(data_dir[i], 
                                if (is.null(sample_name)) NULL else sample_name[i],
                                all_counts[[i]])
        }, BPPARAM = BPPARAM)
      }, error = function(e) {
        warning(paste("Parallel processing failed:", e$message, 
                      "\nFalling back to serial processing"))
        # Fall back to serial processing with preloaded data
        lapply(seq_along(data_dir), function(i) {
          process_one_directory(data_dir[i], 
                                if (is.null(sample_name)) NULL else sample_name[i],
                                all_counts[[i]])
        })
      })
    }
  } else {
    # Serial processing
    objs <- lapply(seq_along(data_dir), function(i) {
      process_one_directory(data_dir[i], 
                            if (is.null(sample_name)) NULL else sample_name[i],
                            NULL)  # No preloaded data for serial processing
    })
  }
  
  # Set names for the list
  if (is.null(sample_name)) {
    names(objs) <- sapply(data_dir, function(x) basename(normalizePath(x)))
  } else {
    names(objs) <- sample_name
  }
  
  message(paste("Successfully created", length(objs), "Seurat objects"))
  
  return(objs)
}

# =============================================================================
# UTILITY FUNCTIONS (kept from original for completeness)
# =============================================================================

# Function to detect data format in a directory
detect_data_format <- function(dir_path) {
  files <- list.files(dir_path, recursive = FALSE)
  
  # Check for Parse files
  parse_terms <- c("cell_metadata", "all_genes", "count_matrix")
  parse_matches <- sapply(parse_terms, function(term) {
    any(grepl(term, files, ignore.case = TRUE))
  })
  
  # Check for 10X files
  tenx_terms <- c("matrix", "features|genes", "barcodes")
  tenx_matches <- sapply(tenx_terms, function(term) {
    any(grepl(term, files, ignore.case = TRUE))
  })
  
  # Check for HDF5 files
  h5_found <- any(grepl("\\.h5(\\.gz|\\.bz2)?$", files, ignore.case = TRUE))
  
  if (sum(parse_matches) >= 2) {
    return("parse")
  } else if (sum(tenx_matches) >= 2 || h5_found) {
    return("tenx")
  } else {
    return("unknown")
  }
}

# Function to find files by pattern
find_file_by_pattern <- function(dir_path, pattern, extensions = c("", ".gz", ".bz2")) {
  files <- list.files(dir_path, full.names = TRUE)
  file_basenames <- basename(files)
  
  for (ext in extensions) {
    search_pattern <- paste0(".*", pattern, ".*", ext, "$")
    matched_files <- files[grepl(search_pattern, file_basenames, ignore.case = TRUE)]
    if (length(matched_files) > 0) {
      return(matched_files[1])
    }
  }
  return(NULL)
}

# Function to decompress file if needed
decompress_file <- function(file_path, temp_dir) {
  if (grepl("\\.gz$", file_path)) {
    decompressed_path <- file.path(temp_dir, gsub("\\.gz$", "", basename(file_path)))
    R.utils::gunzip(file_path, destname = decompressed_path, remove = FALSE, overwrite = TRUE)
    return(decompressed_path)
  } else if (grepl("\\.bz2$", file_path)) {
    decompressed_path <- file.path(temp_dir, gsub("\\.bz2$", "", basename(file_path)))
    R.utils::bunzip2(file_path, destname = decompressed_path, remove = FALSE, overwrite = TRUE)
    return(decompressed_path)
  } else {
    # Copy uncompressed file to temp directory
    temp_path <- file.path(temp_dir, basename(file_path))
    file.copy(file_path, temp_path, overwrite = TRUE)
    return(temp_path)
  }
}

# Function to create and manage temporary directory
create_temp_directory <- function(sample_name) {
  temp_dir <- file.path(tempdir(), paste0(sample_name, "_temp"))
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Return both directory and cleanup function
  cleanup <- function() {
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
    }
  }
  
  return(list(dir = temp_dir, cleanup = cleanup))
}

# Function to validate sample name
get_sample_name <- function(data_dir, sample_name = NULL) {
  if (is.null(sample_name)) {
    sample_name <- basename(dirname(normalizePath(data_dir)))
  }
  return(sample_name)
}

# Function to load Parse Biosciences data
load_parse_data <- function(data_dir, temp_dir, ...) {
  message("Loading Parse Biosciences data...")
  
  # Find required files
  cell_metadata_file <- find_file_by_pattern(data_dir, "cell_metadata")
  gene_metadata_file <- find_file_by_pattern(data_dir, "all_genes")
  count_matrix_file <- find_file_by_pattern(data_dir, "count_matrix")
  
  if (is.null(cell_metadata_file) || is.null(gene_metadata_file) || is.null(count_matrix_file)) {
    stop("Could not find required Parse files (cell_metadata, all_genes, count_matrix)")
  }
  
  # Prepare files
  decompress_file(cell_metadata_file, temp_dir)
  decompress_file(gene_metadata_file, temp_dir)
  decompress_file(count_matrix_file, temp_dir)
  
  # Rename files to standard Parse names
  file.rename(
    file.path(temp_dir, gsub("\\.(gz|bz2)$", "", basename(cell_metadata_file))),
    file.path(temp_dir, "cell_metadata.csv")
  )
  file.rename(
    file.path(temp_dir, gsub("\\.(gz|bz2)$", "", basename(gene_metadata_file))),
    file.path(temp_dir, "all_genes.csv")
  )
  file.rename(
    file.path(temp_dir, gsub("\\.(gz|bz2)$", "", basename(count_matrix_file))),
    file.path(temp_dir, "count_matrix.mtx")
  )
  
  # Load data - need to load Seurat for ReadParseBio
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for loading Parse data")
  }
  counts_data <- Seurat::ReadParseBio(data.dir = temp_dir, ...)
  return(counts_data)
}

# Function to load 10X HDF5 data
load_10x_h5_data <- function(h5_file, temp_dir) {
  message("Found HDF5 file, loading...")
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for loading 10X H5 data")
  }
  
  if (grepl("\\.(gz|bz2)$", h5_file)) {
    # Decompress HDF5 file
    decompressed_h5 <- decompress_file(h5_file, temp_dir)
    counts_data <- Seurat::Read10X_h5(decompressed_h5)
  } else {
    counts_data <- Seurat::Read10X_h5(h5_file)
  }
  
  return(counts_data)
}

# Function to load 10X MTX data
load_10x_mtx_data <- function(data_dir, temp_dir) {
  message("Loading MTX format...")
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for loading 10X MTX data")
  }
  
  # Find required files
  matrix_file <- find_file_by_pattern(data_dir, "matrix")
  features_file <- find_file_by_pattern(data_dir, "features|genes")
  barcodes_file <- find_file_by_pattern(data_dir, "barcodes")
  
  if (is.null(matrix_file) || is.null(features_file) || is.null(barcodes_file)) {
    stop("Could not find required 10X files (matrix, features/genes, barcodes)")
  }
  
  # Copy files to temp directory with standard names
  file.copy(matrix_file, file.path(temp_dir, "matrix.mtx.gz"), overwrite = TRUE)
  file.copy(features_file, file.path(temp_dir, "features.tsv.gz"), overwrite = TRUE)
  file.copy(barcodes_file, file.path(temp_dir, "barcodes.tsv.gz"), overwrite = TRUE)
  
  # Load data
  counts_data <- Seurat::Read10X(data.dir = temp_dir)
  return(counts_data)
}

# Function to load 10X data (handles both HDF5 and MTX)
load_10x_data <- function(data_dir, temp_dir) {
  message("Loading 10X data...")
  
  # Check for HDF5 files first
  h5_file <- find_file_by_pattern(data_dir, "\\.h5")
  
  if (!is.null(h5_file)) {
    return(load_10x_h5_data(h5_file, temp_dir))
  } else {
    return(load_10x_mtx_data(data_dir, temp_dir))
  }
}

# Main data loading function
get_count_matrix <- function(data_dir, data_type = "auto", sample_name = NULL, ...) {
  
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  sample_name <- get_sample_name(data_dir, sample_name)
  
  message(paste("Loading data from:", data_dir))
  
  # Detect format if auto
  if (data_type == "auto") {
    data_type <- detect_data_format(data_dir)
    message(paste("Auto-detected format:", data_type))
  }
  
  # Create temporary directory
  temp_setup <- create_temp_directory(sample_name)
  temp_dir <- temp_setup$dir
  cleanup <- temp_setup$cleanup
  
  # Ensure cleanup happens even if function fails
  on.exit(cleanup())
  
  # Load data based on type
  if (data_type == "parse") {
    counts_data <- load_parse_data(data_dir, temp_dir, ...)
  } else if (data_type == "tenx") {
    counts_data <- load_10x_data(data_dir, temp_dir)
  } else {
    stop(paste("Unsupported or unknown data format:", data_type))
  }
  
  message(paste("Successfully loaded", ncol(counts_data), "cells and", nrow(counts_data), "genes"))
  return(counts_data)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Single directory
# seurat_obj <- dir2seurat("/path/to/sample1/")

# Multiple directories
# data_dirs <- c("/path/to/sample1/", "/path/to/sample2/", "/path/to/sample3/")
# seurat_list <- dir2seurat(data_dirs, sample_name = c("Sample1", "Sample2", "Sample3"))

# With specific parameters
# seurat_obj <- dir2seurat("/path/to/sample1/", 
#                          data_type = "tenx",
#                          min_cells = 5,
#                          min_features = 300,
#                          project = "MyExperiment")