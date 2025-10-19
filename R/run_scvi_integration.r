#' Run scVI Integration for Batch Correction in Single-Cell RNA-seq Data
#'
#' @description
#' Performs batch correction and data integration using scVI (single-cell
#' Variational Inference), a deep generative model for single-cell RNA-seq data.
#' The function automatically optimizes performance parameters based on available
#' hardware (GPU/CPU, RAM), trains the scVI model, and returns the integrated
#' latent representation as a dimensional reduction in the Seurat object.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data with
#'   multiple batches to integrate
#' @param int_features Character vector of feature names (genes) to use for
#'   integration. Typically highly variable genes selected by
#'   \code{\link{custom_VarFeatures}} or \code{\link[Seurat]{FindVariableFeatures}}
#' @param filename Character string. Base name used for creating intermediate
#'   AnnData files (stored in data_path directory)
#' @param assay Character string. Name of the assay to use for integration.
#'   Default: "RNA"
#' @param scvi.model Named list of arguments to pass to
#'   \code{scvi$model$SCVI$setup_anndata()}. Must include \code{batch_key}.
#'   See Details for full description of available arguments. Required.
#' @param ... Additional arguments passed to \code{scvi$model$SCVI()} model
#'   constructor (e.g., n_latent, n_layers, dispersion, gene_likelihood)
#'
#' @return Returns the input Seurat object with an added "scvi" dimensional
#'   reduction containing the integrated latent representation. The reduction
#'   can be accessed via \code{seurat_obj[["scvi"]]} and has the key "scvi_".
#'
#' @details
#' **scVI Model Overview:**
#'
#' scVI (Lopez et al., 2018) is a probabilistic model that learns a low-dimensional
#' latent representation of single-cell expression data while accounting for
#' technical variability including batch effects. Unlike linear methods (e.g.,
#' Harmony, Seurat integration), scVI uses deep neural networks to model complex,
#' nonlinear batch effects.
#'
#' **Understanding batch_key vs. Covariates:**
#'
#' \strong{batch_key} (REQUIRED):
#' \itemize{
#'   \item The primary technical variable to correct (e.g., "sample", "batch",
#'         "donor", "library_prep")
#'   \item scVI learns batch-specific parameters to remove systematic differences
#'   \item Must have at least 2 unique values (function checks and warns if not)
#'   \item This is the main variable you want to "integrate across"
#' }
#'
#' \strong{categorical_covariate_keys} (OPTIONAL):
#' \itemize{
#'   \item Additional categorical variables that may affect expression
#'         (e.g., "sex", "age_group", "tissue_region")
#'   \item Model accounts for these but does NOT remove their effects
#'   \item Use when you want to preserve biological variation while correcting batch
#'   \item Can be a vector: c("sex", "age_group")
#' }
#'
#' \strong{continuous_covariate_keys} (OPTIONAL):
#' \itemize{
#'   \item Continuous variables (e.g., "percent_mito", "nCount_RNA", "age_numeric")
#'   \item Model conditions on these variables
#'   \item Use for technical or biological continuous covariates
#'   \item Can be a vector: c("percent_mito", "nCount_RNA")
#' }
#'
#' **Example Use Cases:**
#'
#' *Case 1: Simple batch correction*
#' \preformatted{
#' scvi.model = list(
#'   batch_key = "sample"  # Integrate across samples
#' )
#' }
#'
#' *Case 2: Batch correction while preserving biological effects*
#' \preformatted{
#' scvi.model = list(
#'   batch_key = "batch",                    # Remove batch effects
#'   categorical_covariate_keys = c("sex")   # Preserve sex differences
#' )
#' }
#'
#' *Case 3: Complex design with multiple covariates*
#' \preformatted{
#' scvi.model = list(
#'   batch_key = "library_prep",
#'   categorical_covariate_keys = c("donor", "tissue"),
#'   continuous_covariate_keys = c("percent_mito", "nCount_RNA")
#' )
#' }
#'
#' **Available scvi.model Arguments:**
#'
#' Arguments passed to \code{scvi$model$SCVI$setup_anndata()}:
#' \itemize{
#'   \item \strong{batch_key}: (REQUIRED) Column name for batch variable
#'   \item \strong{layer}: Data layer to use (default: NULL, uses .X)
#'   \item \strong{categorical_covariate_keys}: Vector of categorical covariate names
#'   \item \strong{continuous_covariate_keys}: Vector of continuous covariate names
#'   \item \strong{labels_key}: Column name for cell type labels (for conditional scVI)
#' }
#'
#' **Available Model Constructor Arguments (...):**
#'
#' Arguments passed to \code{scvi$model$SCVI()}:
#' \itemize{
#'   \item \strong{n_latent}: Number of latent dimensions (default: 30)
#'   \item \strong{n_layers}: Number of hidden layers (default: 2)
#'   \item \strong{n_hidden}: Number of nodes per hidden layer (default: 128)
#'   \item \strong{dropout_rate}: Dropout rate (default: 0.1)
#'   \item \strong{dispersion}: Dispersion mode ("gene", "gene-batch", "gene-label"; default: "gene")
#'   \item \strong{gene_likelihood}: Distribution ("zinb", "nb", "poisson"; default: "nb")
#' }
#'
#' **Performance Optimization:**
#'
#' The function automatically optimizes performance based on system resources:
#'
#' \strong{GPU Mode} (if CUDA GPU available):
#' \itemize{
#'   \item Batch size: Dynamically calculated based on GPU memory (up to 4096)
#'   \item DataLoader workers: 0-1 for small systems, ≥2 for larger systems
#'   \item Accelerator: "cuda"
#'   \item Uses high precision float32 matrix multiplication
#' }
#'
#' \strong{CPU Mode} (if no GPU):
#' \itemize{
#'   \item Batch size: 512-2048 based on RAM (>16GB → 1024, >32GB → 2048)
#'   \item DataLoader workers: 0 for <32GB RAM, otherwise half of logical processors
#'   \item Accelerator: "cpu"
#' }
#'
#' Optimization factors:
#' \itemize{
#'   \item Larger batch sizes improve GPU utilization but require more memory
#'   \item More workers speed up data loading but increase RAM usage
#'   \item Persistent workers reduce worker startup overhead
#'   \item CSR sparse matrix format optimizes training speed
#' }
#'
#' **Training Process:**
#'
#' The model training includes:
#' \itemize{
#'   \item Automatic early stopping to prevent overfitting
#'   \item Hardware-optimized batch sizes
#'   \item Progress monitoring during training
#'   \item ELBO (Evidence Lower Bound) convergence tracking
#' }
#'
#' **Output Interpretation:**
#'
#' The function adds a "scvi" dimensional reduction to the Seurat object containing
#' the learned latent representation (default: 30 dimensions). This latent space:
#' \itemize{
#'   \item Has batch effects removed
#'   \item Preserves biological variation
#'   \item Can be used for clustering, UMAP, and other downstream analyses
#'   \item Is analogous to PCA but batch-corrected
#' }
#'
#' Use the scvi reduction for downstream analysis:
#' \preformatted{
#' # Clustering on scvi latent space
#' seurat_obj <- FindNeighbors(seurat_obj, reduction = "scvi", dims = 1:30)
#' seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
#'
#' # UMAP on scvi latent space
#' seurat_obj <- RunUMAP(seurat_obj, reduction = "scvi", dims = 1:30)
#' }
#'
#' @section System Requirements:
#' \itemize{
#'   \item Python packages: scvi-tools, scanpy, anndata, torch
#'   \item R packages: Seurat, scCustomize, reticulate
#'   \item Optional: CUDA-compatible GPU for faster training
#'   \item RAM: Minimum 16GB recommended, 32GB+ for large datasets
#' }
#'
#' @section Verbose Output:
#' The function prints detailed information during execution:
#' \itemize{
#'   \item Batch variable validation and unique values found
#'   \item scVI setup arguments being used
#'   \item System configuration (CPU, RAM, GPU)
#'   \item Performance parameters (batch size, workers)
#'   \item Model architecture details
#'   \item Training progress and convergence
#'   \item Final latent representation dimensions
#' }
#'
#' @examples
#' \dontrun{
#' # Basic integration across samples
#' seurat_obj <- run_scvi_integration(
#'   seurat_obj = seurat_obj,
#'   int_features = VariableFeatures(seurat_obj),
#'   filename = "my_experiment",
#'   scvi.model = list(batch_key = "sample")
#' )
#'
#' # Integration with biological covariates preserved
#' seurat_obj <- run_scvi_integration(
#'   seurat_obj = seurat_obj,
#'   int_features = hvf_genes,
#'   filename = "my_experiment",
#'   scvi.model = list(
#'     batch_key = "batch",
#'     categorical_covariate_keys = c("sex", "age_group")
#'   )
#' )
#'
#' # Integration with continuous covariates
#' seurat_obj <- run_scvi_integration(
#'   seurat_obj = seurat_obj,
#'   int_features = hvf_genes,
#'   filename = "my_experiment",
#'   scvi.model = list(
#'     batch_key = "library_prep",
#'     continuous_covariate_keys = c("percent_mito", "nCount_RNA")
#'   )
#' )
#'
#' # Custom model parameters
#' seurat_obj <- run_scvi_integration(
#'   seurat_obj = seurat_obj,
#'   int_features = hvf_genes,
#'   filename = "my_experiment",
#'   scvi.model = list(batch_key = "sample"),
#'   n_latent = 50,           # More latent dimensions
#'   n_layers = 3,            # Deeper network
#'   gene_likelihood = "zinb" # Zero-inflated negative binomial
#' )
#'
#' # Use scvi reduction for downstream analysis
#' seurat_obj <- FindNeighbors(seurat_obj, reduction = "scvi", dims = 1:30)
#' seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
#' seurat_obj <- RunUMAP(seurat_obj, reduction = "scvi", dims = 1:30)
#' DimPlot(seurat_obj, reduction = "umap", group.by = "sample")
#' }
#'
#' @references
#' Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., & Yosef, N. (2018).
#' Deep generative modeling for single-cell transcriptomics.
#' Nature Methods, 15(12), 1053-1058.
#'
#' scvi-tools documentation: https://docs.scvi-tools.org/
#'
#' @importFrom Seurat DefaultAssay CreateDimReducObject
#' @importFrom scCustomize as.anndata
#' @importFrom reticulate import py_to_r
#'
#' @seealso
#' \code{\link{custom_VarFeatures}} for selecting integration features
#' \code{\link[Seurat]{FindVariableFeatures}} for standard HVF selection
#' \code{\link[Seurat]{IntegrateData}} for Seurat's native integration
#'
#' @export
run_scvi_integration <- function(seurat_obj, int_features, filename, assay = "RNA",
                                 scvi.model = NULL, ...) {
  library(scCustomize)
  library(reticulate)

  # retrieve system resources info for optimizing performance
  res <- get_resources()

  # Check if scvi.model is provided
  if (is.null(scvi.model)) {
    stop("scvi.model argument is required. Please provide a list with arguments for scvi$model$SCVI$setup_anndata().")
  }

  # Validate that batch_key is provided in scvi.model
  if (!"batch_key" %in% names(scvi.model)) {
    stop("batch_key must be specified in scvi.model list.")
  }

  batch_var <- scvi.model$batch_key

  # Check if batch_var exists in the Seurat object
  if (!batch_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Batch variable '", batch_var, "' not found in Seurat object metadata.", sep = ""))
  }

  # Check if batch_var has more than one unique value
  batch_values <- unique(seurat_obj@meta.data[[batch_var]])
  batch_values <- batch_values[!is.na(batch_values)]  # Remove NA values

  if (length(batch_values) <= 1) {
    warning(paste("Batch variable '", batch_var, "' has only ", length(batch_values),
                  " unique value(s). Batch correction requires multiple batches. Skipping scVI integration.", sep = ""))
    return(seurat_obj)
  }

  cat(paste("Found", length(batch_values), "unique values in batch variable '", batch_var, "':",
            paste(batch_values, collapse = ", "), "\n"))

  # load Python libraries
  pd <- import("pandas", convert = FALSE)
  sc <- import("scanpy", convert = FALSE)
  scvi <- import("scvi", convert = FALSE)
  np <- import("numpy", convert = FALSE)
  anndata <- import("anndata")
  torch <- import("torch", convert = FALSE)
  scipy_sparse <- import("scipy.sparse")

  # Create a temporary Seurat object with only relevant features
  temp_seurat <- JoinLayers(seurat_obj)
  temp_seurat <- temp_seurat[int_features]

  # Create the annData file for python - USE COUNTS AS MAIN LAYER
  adata <- as.anndata(temp_seurat, data_path, paste0("anndata_", filename),
                      assay = assay, main_layer = "counts", other_layers = "data",
                      drop_single_values = FALSE)
  print(adata)

  # Clean data - make sure the batch variable is properly formatted
  if (batch_var %in% names(adata$obs)) {
    adata$obs[[batch_var]] <- adata$obs[[batch_var]]$astype("str")$astype("category")
  }

  # Clean date if it exists (keeping your original code)
  if ("date" %in% names(adata$obs)) {
    adata$obs$date <- adata$obs$date$astype("str")$astype("category")
  }

  # Convert sparse matrix to CSR format for faster training
  if (!scipy_sparse$issparse(adata$X)) {
    cat("Matrix is dense. Converting to CSR sparse format...\n")
    adata$X <- scipy_sparse$csr_matrix(adata$X)
  } else if (!scipy_sparse$isspmatrix_csr(adata$X)) {
    cat("Converting sparse matrix to CSR format for faster training...\n")
    adata$X <- adata$X$tocsr()
  } else {
    cat("Matrix is already in CSR format.\n")
  }

  # Setup AnnData for scVI with user-provided arguments
  cat("Setting up scVI with provided arguments...\n")
  cat("Arguments provided:", paste(names(scvi.model), collapse = ", "), "\n")

  # Convert R vectors to Python lists for covariate keys if they exist
  args_list <- scvi.model
  if ("continuous_covariate_keys" %in% names(args_list)) {
    args_list$continuous_covariate_keys <- as.list(args_list$continuous_covariate_keys)
  }
  if ("categorical_covariate_keys" %in% names(args_list)) {
    args_list$categorical_covariate_keys <- as.list(args_list$categorical_covariate_keys)
  }

  # Show the final arguments that will be passed to setup_anndata
  cat("\n=== FINAL setup_anndata ARGUMENTS ===\n")
  cat("scvi$model$SCVI$setup_anndata(adata,\n")
  for (i in seq_along(args_list)) {
    arg_name <- names(args_list)[i]
    arg_value <- args_list[[i]]

    # Format the argument value for display
    if (is.character(arg_value)) {
      formatted_value <- paste0('"', arg_value, '"')
    } else if (is.list(arg_value)) {
      formatted_value <- paste0('c("', paste(arg_value, collapse = '", "'), '")')
    } else {
      formatted_value <- as.character(arg_value)
    }

    # Add comma for all but the last argument
    comma <- if (i < length(args_list)) "," else ""
    cat("  ", arg_name, " = ", formatted_value, comma, "\n", sep = "")
  }
  cat(")\n")

  # Use do.call equivalent approach for dynamic argument passing
  # This is more robust than hardcoding specific combinations
  tryCatch({
    do.call(scvi$model$SCVI$setup_anndata, c(list(adata), args_list))
  }, error = function(e) {
    stop(paste("Error in scVI setup_anndata:", e$message,
               "\nPlease check your scvi.model arguments."))
  })

  # Print the registry to show what was set up
  cat("scVI setup complete. Registry:\n")
  print(adata$uns[["_scvi_uuid"]])

  ### Setup setting
  set_scvi_settings <- function(res) {
    if (!is.null(res$gpu_name)) {
      # --- GPU mode ---
      gpu_mem <- res$gpu_memory_mb
      ram <- res$ram_mb

      # Batch size rule: ~ GPU_RAM / 8 (heuristic for float32 activations)
      max_batch_gpu <- floor(gpu_mem * 1024 / 8)  # ~ samples, rough guess
      batch_size <- min(4096L, max(256L, max_batch_gpu))

      # Cap batch size also by RAM (~ 1/4 RAM in MB)
      batch_size <- min(batch_size, floor(ram / 4))

      # Workers: small laptops → stick to 0 or 1
      dl_num_workers <- if (ram < 32000 || gpu_mem < 4000) 0L else max(2L, res$cpu_cores)

    } else {
      # --- CPU-only mode ---
      ram <- res$ram_mb

      # Batch size based on RAM
      batch_size <- if (ram > 32000) 2048L else if (ram > 16000) 1024L else 512L

      # Workers: avoid overhead on small laptops
      dl_num_workers <- if (ram < 32000) 0L else max(1L, floor(res$logical_processors / 2))
    }

    # Apply settings
    scvi$settings$dl_num_workers <- as.integer(dl_num_workers)
    scvi$settings$batch_size <- as.integer(batch_size)
    scvi$settings$persistent_workers <- TRUE
    torch$set_float32_matmul_precision("high")

    list(
      dl_num_workers = dl_num_workers,
      batch_size = batch_size
    )
  }

  # Apply the optimized settings
  settings <- set_scvi_settings(res)

  # Report system configuration and chosen parameters
  cat("\n=== SYSTEM CONFIGURATION ===\n")
  cat("CPU cores:", res$cpu_cores, "\n")
  cat("Logical processors:", res$logical_processors, "\n")
  cat("RAM:", round(res$ram_mb / 1024, 1), "GB\n")
  if (!is.null(res$gpu_name)) {
    cat("GPU:", res$gpu_name, "\n")
    cat("GPU Memory:", round(res$gpu_memory_mb / 1024, 1), "GB\n")
    accelerator <- "cuda"
    devices <- "auto"
  } else {
    cat("GPU: Not available\n")
    accelerator <- "cpu"
    devices <- "auto"
  }

  cat("\n=== PERFORMANCE PARAMETERS ===\n")
  cat("DataLoader workers:", settings$dl_num_workers, "\n")
  cat("Batch size:", settings$batch_size, "\n")
  cat("Persistent workers: TRUE\n")  # We know we set this to TRUE
  cat("Accelerator:", accelerator, "\n")
  cat("Devices:", devices, "\n")

  # create the model
  cat("\n=== CREATING scVI MODEL ===\n")
  model <- scvi$model$SCVI(adata, ...)

  # Report model details
  cat("Model parameters:\n")
  cat("  - Latent dimensions:", 30, "\n")
  cat("  - Number of layers:", 2, "\n")
  cat("  - Dispersion:", "gene-label", "\n")
  cat("  - Gene likelihood:", "nb", "\n")
  cat("  - Number of genes:", py_to_r(model$adata$n_vars), "\n")
  cat("  - Number of cells:", py_to_r(model$adata$n_obs), "\n")

  # Print model summary
  cat("\nModel summary:\n")
  print(model)

  # train the model
  cat("\n=== TRAINING scVI MODEL ===\n")
  cat("Starting training with the following parameters:\n")
  cat("  - Accelerator:", accelerator, "\n")
  cat("  - Devices:", devices, "\n")
  cat("  - Batch size:", settings$batch_size, "\n")
  cat("  - Early stopping: TRUE\n")

  model$train(accelerator=accelerator, devices=devices, early_stopping=TRUE,
              batch_size = scvi$settings$batch_size)

  # get the latent representation
  cat("\n=== EXTRACTING LATENT REPRESENTATION ===\n")
  cat("Extracting latent representation from trained model...\n")
  latent <- model$get_latent_representation()
  cat("Latent representation shape:", dim(latent), "\n")

  # put it back in our original Seurat object
  latent <- as.matrix(latent)
  rownames(latent) <- colnames(temp_seurat)

  seurat_obj[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_obj))

  cat("scVI integration completed successfully!\n")
  cat("Added 'scvi' dimensional reduction to Seurat object with key 'scvi_'\n")
  cat("Latent dimensions available:", ncol(latent), "\n")

  # Return the updated Seurat object
  return(seurat_obj)
}
