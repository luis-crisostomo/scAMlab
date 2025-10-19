# List available metadata
# list_metadata <- function(seurat_obj) {
#   if (!inherits(seurat_obj, "Seurat")) {
#     stop("Input must be a Seurat object")
#   }
#   
#   available_cols <- colnames(seurat_obj@meta.data)
#   cat("Available metadata columns for filtering:\n")
#   for (i in seq_along(available_cols)) {
#     col_name <- available_cols[i]
#     col_class <- class(seurat_obj@meta.data[[col_name]])[1]
#     col_summary <- if (is.numeric(seurat_obj@meta.data[[col_name]])) {
#       paste0("(", col_class, ", range: ", 
#              round(min(seurat_obj@meta.data[[col_name]], na.rm = TRUE), 2), "-",
#              round(max(seurat_obj@meta.data[[col_name]], na.rm = TRUE), 2), ")")
#     } else {
#       paste0("(", col_class, ")")
#     }
#     cat(paste0(i, ". ", col_name, " ", col_summary, "\n"))
#   }
#   return(available_cols)
# }

#' @export
list_metadata <- function(seurat_obj, show_common_only = FALSE) {
  # Helper function to process a single Seurat object
  process_single_seurat <- function(obj, obj_name = NULL) {
    if (!inherits(obj, "Seurat")) {
      stop(paste("Object", ifelse(is.null(obj_name), "", paste0("'", obj_name, "'")), 
                 "is not a Seurat object"))
    }
    
    available_cols <- colnames(obj@meta.data)
    
    # Create summary for each column
    col_info <- lapply(available_cols, function(col_name) {
      col_data <- obj@meta.data[[col_name]]
      col_class <- class(col_data)[1]
      
      if (is.numeric(col_data)) {
        col_summary <- paste0("(", col_class, ", range: ", 
                              round(min(col_data, na.rm = TRUE), 2), "-",
                              round(max(col_data, na.rm = TRUE), 2), ")")
      } else if (is.factor(col_data) || is.character(col_data)) {
        unique_vals <- length(unique(col_data))
        col_summary <- paste0("(", col_class, ", ", unique_vals, " unique values)")
      } else {
        col_summary <- paste0("(", col_class, ")")
      }
      
      list(name = col_name, class = col_class, summary = col_summary)
    })
    
    return(list(columns = available_cols, info = col_info))
  }
  
  # Check if input is a list
  if (is.list(seurat_obj) && !inherits(seurat_obj, "Seurat")) {
    # Handle list of Seurat objects
    if (length(seurat_obj) == 0) {
      stop("Input list is empty")
    }
    
    # Get names for the list elements
    obj_names <- names(seurat_obj)
    if (is.null(obj_names)) {
      obj_names <- paste0("Object_", seq_along(seurat_obj))
    }
    
    # Process each object
    all_results <- list()
    all_columns <- list()
    
    for (i in seq_along(seurat_obj)) {
      cat("\n=== ", obj_names[i], " ===\n")
      result <- process_single_seurat(seurat_obj[[i]], obj_names[i])
      all_results[[obj_names[i]]] <- result
      all_columns[[i]] <- result$columns
      
      # Display results for this object
      cat("Available metadata columns:\n")
      for (j in seq_along(result$info)) {
        info <- result$info[[j]]
        cat(paste0(j, ". ", info$name, " ", info$summary, "\n"))
      }
    }
    
    # Show common columns across all objects if requested
    if (show_common_only && length(seurat_obj) > 1) {
      common_cols <- Reduce(intersect, all_columns)
      cat("\n=== COMMON COLUMNS ACROSS ALL OBJECTS ===\n")
      if (length(common_cols) > 0) {
        for (i in seq_along(common_cols)) {
          cat(paste0(i, ". ", common_cols[i], "\n"))
        }
      } else {
        cat("No common columns found across all objects.\n")
      }
      return(invisible(list(all_metadata = all_results, common_columns = common_cols)))
    }
    
    # Return just the column names for each object (similar to single object behavior)
    result_cols <- lapply(all_results, function(x) x$columns)
    return(invisible(result_cols))
    
  } else {
    # Handle single Seurat object
    result <- process_single_seurat(seurat_obj)
    
    cat("Available metadata columns for filtering:\n")
    for (i in seq_along(result$info)) {
      info <- result$info[[i]]
      cat(paste0(i, ". ", info$name, " ", info$summary, "\n"))
    }
    
    return(invisible(result$columns))
  }
}