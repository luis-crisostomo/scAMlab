#' Calculate Comprehensive QC Metrics for Single-Cell RNA-seq Data
#'
#' @description
#' Computes quality control metrics combining standard Seurat tutorial metrics
#' with additional QC measures recommended by the Babraham Institute. Supports
#' both single Seurat objects and lists of Seurat objects with optional parallel
#' processing and memory management.
#'
#' @param seurat_obj A Seurat object or a list of Seurat objects
#' @param assay Character string. Name of the assay to use for QC calculations.
#'   Default: "RNA"
#' @param max_memory_gb Numeric or NULL. Maximum memory usage threshold in GB.
#'   If specified, the function will stop if memory usage exceeds this limit.
#'   Default: NULL (no limit)
#' @param chunk_size Integer. Number of cells to process in each chunk for
#'   memory-efficient computation of the largest gene per cell. Default: 1000
#' @param exclude_genes Character vector or NULL. Gene names or patterns to
#'   exclude from largest gene calculations (e.g., mitochondrial or ribosomal
#'   genes). Default: NULL (no exclusions)
#' @param mt.fraction Character string. Either a regex pattern to identify
#'   mitochondrial genes (e.g., "^mt-|^MT-") or a character vector of
#'   mitochondrial gene names. Default: "^mt-|^MT-"
#' @param rb.fraction Character string. Either a regex pattern to identify
#'   ribosomal genes (e.g., "^Rp\[sl\]|^RP\[SL\]") or a character vector of
#'   ribosomal gene names. Default: "^Rp\[sl\]|^RP\[SL\]"
#' @param BPPARAM BiocParallelParam object for parallel processing when
#'   computing the largest gene per cell (e.g., MulticoreParam, SnowParam).
#'   If NULL, uses SerialParam for sequential processing. Default: NULL
#'
#' @return
#' \itemize{
#'   \item For single Seurat object: Returns the Seurat object with added QC
#'         metrics in the metadata
#'   \item For list of Seurat objects: Returns a named list of Seurat objects,
#'         each with added QC metrics
#' }
#'
#' @details
#' This function calculates QC metrics from two established workflows:
#'
#' **Standard Seurat Tutorial Metrics**
#' (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html):
#' \itemize{
#'   \item nCount_RNA: Total UMI counts per cell
#'   \item nFeature_RNA: Number of genes detected per cell
#'   \item percent.mt: Percentage of mitochondrial gene expression
#' }
#'
#' **Babraham Institute Additional Metrics**
#' (https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html):
#' \itemize{
#'   \item log10GenesPerUMI: Novelty score (log10(nFeature) / log10(nCount))
#'   \item percent.rb: Percentage of ribosomal gene expression
#'   \item largest_gene: Gene with highest expression in each cell
#'   \item percent.Largest.Gene: Percentage of counts from the largest gene
#'   \item novelty_diff: Deviation from expected log10(genes) based on
#'         log10(UMIs), calculated as residuals from linear regression
#' }
#'
#' @section Gene Pattern Matching:
#' The function accepts flexible input for mitochondrial and ribosomal genes:
#' \itemize{
#'   \item **Regex patterns** (detected by presence of regex special characters):
#'         Uses pattern matching to identify genes
#'   \item **Gene name vectors**: Directly uses provided gene names
#' }
#'
#' Common patterns:
#' \itemize{
#'   \item Mouse mitochondrial: "^mt-" or "^Mt-"
#'   \item Human mitochondrial: "^MT-"
#'   \item Ribosomal proteins: "^Rp\[sl\]" (mouse) or "^RP\[SL\]" (human)
#' }
#'
#' @section Memory Management:
#' The function implements several memory optimization strategies:
#' \itemize{
#'   \item **Chunked processing**: Cells are processed in batches to manage memory
#'   \item **Memory monitoring**: Optional memory limit enforcement
#'   \item **Garbage collection**: Periodic cleanup during computation
#'   \item **Parallel processing**: Distributes memory load across workers
#' }
#'
#' @section Novelty Score:
#' The novelty_diff metric identifies cells with unusual gene detection rates:
#' \itemize{
#'   \item Fits a linear model: log10(nFeature) ~ log10(nCount)
#'   \item Calculates residuals as novelty_diff
#'   \item Positive values: more genes than expected (good complexity)
#'   \item Negative values: fewer genes than expected (potential low-quality cells)
#'   \item The fitted linear model is stored in seurat_obj@tools$novelty.lm
#' }
#'
#' @section Parallel Processing:
#' When calculating the largest gene per cell across many cells, parallel
#' processing can significantly speed up computation:
#' \preformatted{
#' library(BiocParallel)
#' # For Unix-like systems (Linux, macOS)
#' param <- MulticoreParam(workers = 4)
#'
#' # For Windows
#' param <- SnowParam(workers = 4, type = "SOCK")
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' seurat_obj <- run_babraham_qc(seurat_obj)
#'
#' # Custom mitochondrial and ribosomal patterns for human data
#' seurat_obj <- run_babraham_qc(
#'   seurat_obj = seurat_obj,
#'   mt.fraction = "^MT-",
#'   rb.fraction = "^RP\[SL\]"
#' )
#'
#' # Exclude mitochondrial and ribosomal genes from largest gene calculation
#' seurat_obj <- run_babraham_qc(
#'   seurat_obj = seurat_obj,
#'   exclude_genes = c("^mt-", "^Mt-", "^Rp\[sl\]", "^RP\[SL\]")
#' )
#'
#' # Process with memory limit and parallel processing
#' library(BiocParallel)
#' param <- MulticoreParam(workers = 4)
#' seurat_obj <- run_babraham_qc(
#'   seurat_obj = seurat_obj,
#'   max_memory_gb = 32,
#'   chunk_size = 500,
#'   BPPARAM = param
#' )
#'
#' # Process multiple samples
#' seurat_list <- list(
#'   Control = seurat_ctrl,
#'   Treatment = seurat_trt
#' )
#'
#' processed_list <- run_babraham_qc(
#'   seurat_obj = seurat_list,
#'   BPPARAM = MulticoreParam(workers = 4)
#' )
#'
#' # Use specific gene lists instead of patterns
#' mt_genes <- c("mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp6")
#' rb_genes <- c("Rpl3", "Rpl4", "Rps3", "Rps4")
#'
#' seurat_obj <- run_babraham_qc(
#'   seurat_obj = seurat_obj,
#'   mt.fraction = mt_genes,
#'   rb.fraction = rb_genes
#' )
#' }
#'
#' @references
#' Seurat PBMC Tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#'
#' Babraham Institute scRNA-seq Course:
#' https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html
#'
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom Matrix colSums
#' @importFrom stringr str_detect
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @importFrom stats lm
#'
#' @seealso
#' \code{\link{plot_QC_classic}} for visualizing basic QC metrics
#' \code{\link{plot_QC_UMAP}} for UMAP-based QC visualizations
#' \code{\link[BiocParallel]{BiocParallelParam}} for parallel processing options
#'
#' @export
run_babraham_qc <- function(seurat_obj,
                            assay = "RNA",
                            max_memory_gb = NULL,
                            chunk_size = 1000,
                            exclude_genes = NULL,
                            mt.fraction = "^mt-|^MT-",
                            rb.fraction = "^Rp[sl]|^RP[SL]", #|^Rn\\d+s",
                            BPPARAM = NULL) {

  # Load required packages
  require(stringr)
  require(Matrix)
  require(BiocParallel)

  # Handle BiocParallel parameter
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
    message("Using SerialParam() for sequential processing")
  } else {
    message(paste("Using", class(BPPARAM)[1], "with", bpnworkers(BPPARAM), "workers"))
  }

  # Check if input is a list of Seurat objects or a single object
  if (is.list(seurat_obj) && !inherits(seurat_obj, "Seurat")) {
    # Input is a list of Seurat objects
    message(paste("Processing", length(seurat_obj), "Seurat objects..."))

    # Process each object in the list
    processed_objects <- vector("list", length(seurat_obj))
    names(processed_objects) <- names(seurat_obj)

    for (i in seq_along(seurat_obj)) {
      obj_name <- if (is.null(names(seurat_obj)[i]) || names(seurat_obj)[i] == "") {
        paste("Object", i)
      } else {
        names(seurat_obj)[i]
      }

      message(paste("\n--- Processing", obj_name, "---"))

      processed_objects[[i]] <- run_babraham_qc_single(
        seurat_obj = seurat_obj[[i]],
        assay = assay,
        max_memory_gb = max_memory_gb,
        chunk_size = chunk_size,
        exclude_genes = exclude_genes,
        mitochondrial = mt.fraction,
        ribosomal = rb.fraction,
        BPPARAM = BPPARAM
      )
    }

    message("\n=== All objects processed successfully ===")
    return(processed_objects)

  } else if (inherits(seurat_obj, "Seurat")) {
    # Input is a single Seurat object
    return(run_babraham_qc_single(
      seurat_obj = seurat_obj,
      assay = assay,
      max_memory_gb = max_memory_gb,
      chunk_size = chunk_size,
      exclude_genes = exclude_genes,
      mitochondrial = mt.fraction,
      ribosomal = rb.fraction,
      BPPARAM = BPPARAM
    ))
  } else {
    stop("Input must be either a Seurat object or a list of Seurat objects")
  }
}

# Internal function to process a single Seurat object
run_babraham_qc_single <- function(seurat_obj,
                                   assay = "RNA",
                                   max_memory_gb = NULL,
                                   chunk_size = 1000,
                                   exclude_genes = NULL,
                                   mitochondrial = "^mt-|^MT-",
                                   ribosomal = "^Rp[sl]|^RP[SL]|^Rn\\d+s",
                                   BPPARAM = SerialParam()) {

  # Memory monitoring function
  check_memory <- function(max_gb = NULL) {
    if (!is.null(max_gb)) {
      current_usage <- as.numeric(object.size(ls(envir = .GlobalEnv))) / 1e9
      available_memory <- tryCatch({
        if (.Platform$OS.type == "windows") {
          memory.limit() / 1024  # Convert MB to GB on Windows
        } else {
          # On Unix systems, try to get available memory
          system("free -g | awk 'NR==2{print $7}'", intern = TRUE) %>% as.numeric()
        }
      }, error = function(e) NULL)

      if (current_usage > max_gb) {
        stop(paste("Memory usage exceeded limit:", round(current_usage, 2), "GB >", max_gb, "GB"))
      }
    }
  }

  # Validate assay exists
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste("Assay", assay, "not found in Seurat object. Available assays:",
               paste(names(seurat_obj@assays), collapse = ", ")))
  }

  message(paste("Running QC analysis on assay:", assay))

  # Get counts matrix from specified assay
  counts_matrix <- seurat_obj@assays[[assay]]@layers$counts
  if (is.null(counts_matrix)) {
    counts_matrix <- seurat_obj@assays[[assay]]@counts
  }

  # Calculate basic metrics if not present
  ncount_col <- paste0("nCount_", assay)
  nfeature_col <- paste0("nFeature_", assay)

  if (!ncount_col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data[[ncount_col]] <- Matrix::colSums(counts_matrix)
  }

  if (!nfeature_col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data[[nfeature_col]] <- Matrix::colSums(counts_matrix > 0)
  }

  # Calculate log10GenesPerUMI
  seurat_obj@meta.data$log10GenesPerUMI <- log10(seurat_obj@meta.data[[nfeature_col]]) /
    log10(seurat_obj@meta.data[[ncount_col]])

  # Flexible mitochondrial gene detection
  if (!is.null(mitochondrial)) {
    if (is.character(mitochondrial) && length(mitochondrial) == 1 &&
        any(grepl("\\^|\\$|\\[|\\]|\\*|\\+|\\?|\\||\\(|\\)", mitochondrial))) {
      # Input appears to be a regex pattern
      message(paste("Using mitochondrial pattern:", mitochondrial))
      seurat_obj@meta.data$percent.mt <- PercentageFeatureSet(seurat_obj,
                                                              pattern = mitochondrial,
                                                              assay = assay)
    } else if (is.character(mitochondrial)) {
      # Input is a vector of gene names
      mt_genes_present <- intersect(mitochondrial, rownames(seurat_obj))
      if (length(mt_genes_present) == 0) {
        warning("No provided mitochondrial genes found in the dataset")
        seurat_obj@meta.data$percent.mt <- 0
      } else {
        message(paste("Using", length(mt_genes_present), "user-provided mitochondrial genes"))
        seurat_obj@meta.data$percent.mt <- PercentageFeatureSet(seurat_obj,
                                                                features = mt_genes_present,
                                                                assay = assay)
      }
    } else {
      warning("Mitochondrial argument must be a character vector (genes) or pattern")
      seurat_obj@meta.data$percent.mt <- 0
    }
  } else {
    warning("No mitochondrial genes provided")
    seurat_obj@meta.data$percent.mt <- 0
  }

  # Flexible ribosomal gene detection
  if (!is.null(ribosomal)) {
    if (is.character(ribosomal) && length(ribosomal) == 1 &&
        any(grepl("\\^|\\$|\\[|\\]|\\*|\\+|\\?|\\||\\(|\\)", ribosomal))) {
      # Input appears to be a regex pattern
      message(paste("Using ribosomal pattern:", ribosomal))
      seurat_obj@meta.data$percent.rb <- PercentageFeatureSet(seurat_obj,
                                                              pattern = ribosomal,
                                                              assay = assay)
    } else if (is.character(ribosomal)) {
      # Input is a vector of gene names
      rb_genes_present <- intersect(ribosomal, rownames(seurat_obj))
      if (length(rb_genes_present) == 0) {
        warning("No provided ribosomal genes found in the dataset")
        seurat_obj@meta.data$percent.rb <- 0
      } else {
        message(paste("Using", length(rb_genes_present), "user-provided ribosomal genes"))
        seurat_obj@meta.data$percent.rb <- PercentageFeatureSet(seurat_obj,
                                                                features = rb_genes_present,
                                                                assay = assay)
      }
    } else {
      warning("Ribosomal argument must be a character vector (genes) or pattern")
      seurat_obj@meta.data$percent.rb <- 0
    }
  } else {
    warning("No ribosomal genes provided")
    seurat_obj@meta.data$percent.rb <- 0
  }

  # Initial memory check
  check_memory(max_memory_gb)

  message("Starting babraham_qc analysis...")
  initial_memory <- gc()

  # Handle gene exclusion
  if (is.null(exclude_genes)) {
    exclude_genes <- character(0)
    message("No genes will be excluded from analysis")
  } else {
    message(paste("Excluding", length(exclude_genes), "genes:",
                  paste(head(exclude_genes, 5), collapse = ", "),
                  if(length(exclude_genes) > 5) "..." else ""))
  }

  # More efficient filtering
  if (length(exclude_genes) > 0) {
    gene_mask <- !str_detect(rownames(seurat_obj), paste(exclude_genes, collapse = "|"))
  } else {
    gene_mask <- rep(TRUE, nrow(seurat_obj))
  }

  # Work with filtered matrix
  filtered_counts <- counts_matrix[gene_mask, , drop = FALSE]

  check_memory(max_memory_gb)
  gc(verbose = FALSE)

  # BiocParallel chunked processing function
  process_chunk_biocparallel <- function(chunk_indices, counts_matrix) {
    # Load required packages in worker
    require(Matrix)

    chunk_matrix <- as.matrix(counts_matrix[, chunk_indices, drop = FALSE])

    # Find max values and indices for this chunk
    chunk_max <- apply(chunk_matrix, 2, max)
    chunk_which_max <- apply(chunk_matrix, 2, which.max)

    return(list(max_vals = chunk_max, max_indices = chunk_which_max))
  }

  # Setup for chunked processing
  n_cells <- ncol(filtered_counts)
  n_chunks <- ceiling(n_cells / chunk_size)

  # Create chunk indices
  chunk_list <- vector("list", n_chunks)
  for (i in 1:n_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_cells)
    chunk_list[[i]] <- start_idx:end_idx
  }

  message(paste("Processing", n_cells, "cells in", n_chunks, "chunks using",
                class(BPPARAM)[1], "with", bpnworkers(BPPARAM), "workers..."))

  # Use BiocParallel for processing
  if (inherits(BPPARAM, "SerialParam") || bpnworkers(BPPARAM) == 1) {
    # Sequential processing
    results <- vector("list", n_chunks)
    for (i in 1:n_chunks) {
      results[[i]] <- process_chunk_biocparallel(chunk_list[[i]], filtered_counts)

      if (i %% 5 == 0) {
        check_memory(max_memory_gb)
        gc(verbose = FALSE)
        message(paste("Processed chunk", i, "of", n_chunks))
      }
    }
  } else {
    # Parallel processing using BiocParallel
    results <- bplapply(chunk_list, function(chunk_idx) {
      process_chunk_biocparallel(chunk_idx, filtered_counts)
    }, BPPARAM = BPPARAM)
  }

  # Combine results
  largest_count <- numeric(n_cells)
  largest_index <- integer(n_cells)

  for (i in 1:n_chunks) {
    chunk_indices <- chunk_list[[i]]
    largest_count[chunk_indices] <- results[[i]]$max_vals
    largest_index[chunk_indices] <- results[[i]]$max_indices
  }

  # Convert indices to gene names using the filtered gene list
  filtered_genes <- rownames(seurat_obj)[gene_mask]
  largest_gene <- filtered_genes[largest_index]

  # Calculate percentages using the correct count column
  nCount_RNA <- seurat_obj@meta.data[[ncount_col]]
  percent_largest <- 100 * largest_count / nCount_RNA

  # Add to metadata
  seurat_obj@meta.data$largest_gene <- largest_gene
  seurat_obj@meta.data$percent.Largest.Gene <- percent_largest

  check_memory(max_memory_gb)
  gc(verbose = FALSE)

  # Novelty score calculation
  nFeature_RNA <- seurat_obj@meta.data[[nfeature_col]]
  novelty.lm <- lm(log10(nFeature_RNA) ~ log10(nCount_RNA))

  # Calculate novelty difference
  seurat_obj@meta.data$novelty_diff <- log10(nFeature_RNA) -
    ((log10(nCount_RNA) * novelty.lm$coefficients[2]) + novelty.lm$coefficients[1])

  # Store the linear model
  if (is.null(seurat_obj@tools)) {
    seurat_obj@tools <- list()
  }
  seurat_obj@tools$novelty.lm <- novelty.lm

  # Final cleanup and memory report
  final_memory <- gc()
  memory_diff <- final_memory[2,2] - initial_memory[2,2]

  message(paste("Analysis complete. Memory change:", round(memory_diff, 2), "MB"))
  message("QC metrics added:")
  message("  - log10GenesPerUMI")
  message("  - percent.mt")
  message("  - percent.rb")
  message("  - largest_gene")
  message("  - percent.Largest.Gene")
  message("  - novelty_diff")

  return(seurat_obj)
}
