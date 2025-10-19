#' @export

custom_cluster_annotations <- function(seurat_list, cluster_annotations) {
  
  # Get sample names from seurat_list
  sample_names <- names(seurat_list)
  if(is.null(sample_names)) {
    sample_names <- paste0("sample", seq_along(seurat_list))
    names(seurat_list) <- sample_names
  }
  
  # Apply annotations to each sample
  for(sample_name in sample_names) {
    if(sample_name %in% names(cluster_annotations)) {
      
      cat("Processing sample:", sample_name, "\n")
      
      # Get current seurat object
      seurat_obj <- seurat_list[[sample_name]]
      
      # Get current active identity as character
      current_idents <- as.character(Idents(seurat_obj))
      
      # Check what clusters exist in the data vs annotations
      existing_clusters <- unique(current_idents)
      annotation_clusters <- names(cluster_annotations[[sample_name]])
      
      cat("  Existing clusters in data:", paste(existing_clusters, collapse = ", "), "\n")
      cat("  Annotation clusters:", paste(annotation_clusters, collapse = ", "), "\n")
      
      # Find missing clusters
      missing_in_data <- setdiff(annotation_clusters, existing_clusters)
      missing_in_annotations <- setdiff(existing_clusters, annotation_clusters)
      
      if(length(missing_in_data) > 0) {
        cat("  ⚠ Clusters in annotations but not in data:", paste(missing_in_data, collapse = ", "), "\n")
      }
      
      if(length(missing_in_annotations) > 0) {
        cat("  ⚠ Clusters in data but not in annotations:", paste(missing_in_annotations, collapse = ", "), "\n")
      }
      
      # Create cell type vector for all cells
      # Initialize with NA for safety
      cell_types <- rep(NA, length(current_idents))
      names(cell_types) <- names(current_idents)
      
      # Map cluster numbers to cell types
      for(cluster_id in names(cluster_annotations[[sample_name]])) {
        if(cluster_id %in% existing_clusters) {
          # Find cells belonging to this cluster
          cluster_cells <- current_idents == cluster_id
          cell_types[cluster_cells] <- cluster_annotations[[sample_name]][cluster_id]
        }
      }
      
      # Check for any unmapped cells
      unmapped_cells <- sum(is.na(cell_types))
      if(unmapped_cells > 0) {
        cat("  ⚠", unmapped_cells, "cells could not be mapped to cell types\n")
        # You might want to assign these to "Unknown" instead of NA
        cell_types[is.na(cell_types)] <- "Unknown"
      }
      
      # Add cell types to metadata using AddMetaData (safer method)
      seurat_list[[sample_name]] <- AddMetaData(seurat_obj, 
                                                metadata = cell_types, 
                                                col.name = "cell_type")
      
      # Print summary
      cell_type_counts <- table(cell_types)
      cat("  ✓ Cell type distribution:\n")
      for(ct in names(cell_type_counts)) {
        cat("    ", ct, ":", cell_type_counts[ct], "cells\n")
      }
      
    } else {
      cat("⚠ No annotations found for", sample_name, "\n")
      cat("  Available annotation samples:", paste(names(cluster_annotations), collapse = ", "), "\n")
    }
    cat("\n")
  }
  
  return(seurat_list)
}