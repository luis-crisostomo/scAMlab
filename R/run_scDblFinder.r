#' @export
run_scDblFinder <- function(seurat_obj, clusters = TRUE, assay = "RNA", ...){
  
  # Handle list input
  if(is.list(seurat_obj) && !inherits(seurat_obj, "Seurat")) {
    return(run_scDblFinder_list(seurat_obj, clusters = clusters, assay = assay, ...))
  }
  
  # Single object processing (your existing code)
  if(!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object or list of Seurat objects")
  }
  
  if(!is.logical(clusters)) {
    stop("clusters parameter must be TRUE or FALSE")
  }
  
  # Validate assay parameter
  if(!is.character(assay) || length(assay) != 1) {
    stop("assay parameter must be a single character string")
  }
  
  # Check if the specified assay exists
  if(!assay %in% names(seurat_obj@assays)) {
    available_assays <- names(seurat_obj@assays)
    stop(paste("Assay '", assay, "' not found in Seurat object. Available assays: ", 
               paste(available_assays, collapse = ", "), sep = ""))
  }
  
  invisible(gc())
  
  message(paste("Using assay:", assay))
  message("Converting to SingleCellExperiment...")
  
  # Convert to SCE using the specified assay
  sce <- as.SingleCellExperiment(seurat_obj, assay = assay)
  invisible(gc())
  
  if(clusters){
    message("Running fastcluster()...")
    sce$clusters <- fastcluster(sce)
    message("Fastcluster complete!")
    invisible(gc())
  } 
  
  message("Running scDblFinder...")
  sce <- scDblFinder(sce, ...)
  message("scDblFinder complete!")
  invisible(gc())
  
  doublet_scores <- sce$scDblFinder.score
  doublet_class <- sce$scDblFinder.class
  
  seurat_obj$doublet_scores <- doublet_scores
  seurat_obj$doublet_class <- doublet_class 
  
  if(clusters && "clusters" %in% colnames(colData(sce))) {
    seurat_obj$fastcluster_clusters <- sce$clusters
  }
  
  rm(sce)
  invisible(gc())
  
  n_doublets <- sum(doublet_class == "doublet", na.rm = TRUE)
  n_cells <- length(doublet_class)
  
  message(paste("Detected", n_doublets, "doublets out of", n_cells, "cells."))
  
  return(seurat_obj)
}

# Updated helper function for list processing
run_scDblFinder_list <- function(seurat_list, clusters = TRUE, assay = "RNA", ...) {
  
  # Validate list input
  if(!all(sapply(seurat_list, function(x) inherits(x, "Seurat")))) {
    stop("All elements in the list must be Seurat objects")
  }
  
  n_objects <- length(seurat_list)
  
  # Get additional arguments
  dots <- list(...)
  
  # Improved helper function to expand parameters
  expand_param <- function(param, n, param_name) {
    if(length(param) == 1) {
      # For single values, replicate appropriately based on type
      if(inherits(param, "BiocParallelParam") || isS4(param)) {
        # For S4 objects, return as list
        return(replicate(n, param, simplify = FALSE))
      } else {
        # For regular objects, use rep
        return(rep(param, n))
      }
    } else if(length(param) == n) {
      # Already correct length
      if(is.list(param)) {
        return(param)
      } else {
        return(as.list(param))
      }
    } else {
      stop(paste("Parameter '", param_name, "' must have length 1 or", n, "but has length", length(param)))
    }
  }
  
  # Expand clusters parameter
  clusters_expanded <- expand_param(clusters, n_objects, "clusters")
  
  # Expand assay parameter
  assay_expanded <- expand_param(assay, n_objects, "assay")
  
  # Expand all additional parameters
  dots_expanded <- mapply(function(param, param_name) {
    expand_param(param, n_objects, param_name)
  }, dots, names(dots), SIMPLIFY = FALSE)
  
  # Process each object
  message(paste("Processing", n_objects, "Seurat objects..."))
  
  results <- vector("list", n_objects)
  names(results) <- names(seurat_list)
  
  for(i in seq_along(seurat_list)) {
    
    message(paste("\n--- Processing object", i, "of", n_objects, 
                  if(!is.null(names(seurat_list)[i])) paste("(", names(seurat_list)[i], ")") else "", "---"))
    
    # Extract parameters for this iteration
    current_dots <- lapply(dots_expanded, function(param) {
      if(is.list(param)) {
        param[[i]]
      } else {
        param[i]
      }
    })
    
    # Get the assay for this iteration
    current_assay <- if(is.list(assay_expanded)) assay_expanded[[i]] else assay_expanded[i]
    
    # Run scDblFinder on current object
    results[[i]] <- do.call(run_scDblFinder, 
                            c(list(seurat_obj = seurat_list[[i]], 
                                   clusters = if(is.list(clusters_expanded)) clusters_expanded[[i]] else clusters_expanded[i],
                                   assay = current_assay), 
                              current_dots))
    
    # Memory cleanup between objects
    invisible(gc())
  }
  
  message("\n=== All objects processed successfully! ===")
  
  return(results)
}