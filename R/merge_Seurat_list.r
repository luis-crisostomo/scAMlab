#' @export
merge_Seurat_list <- function(seurat_list,
                              project = getOption(x = "Seurat.object.project", default = "SeuratProject"),
                              join_layers = FALSE,  # Changed default to FALSE
                              calculate_inflections = TRUE) {
  
  # Check if list is named
  if (is.null(names(seurat_list))) {
    names(seurat_list) <- paste0("Sample_", seq_along(seurat_list))
    warning("Seurat list was not named. Using default names: Sample_1, Sample_2, etc.")
  }
  
  message("Merging ", length(seurat_list), " Seurat objects...")
  
  # Merge objects while preserving sample identity
  mergedSeurat <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),  # Sample prefixes
    merge.data = TRUE,
    project = project
  )
  
  # Add sample information to metadata for batch effect analysis
  message("Adding sample metadata...")
  sample_ids <- rep(names(seurat_list), sapply(seurat_list, ncol))
  mergedSeurat$sample_id <- sample_ids
  
  # Optional: Join layers only if specifically requested
  # (Better to keep separate for integration methods)
  if (join_layers) {
    message("Joining layers...")
    mergedSeurat <- JoinLayers(mergedSeurat)
  } else {
    message("Keeping layers separate for integration analysis")
    message("Available layers: ", paste(Layers(mergedSeurat), collapse = ", "))
  }
  
  # Calculate barcode inflections if requested
  if (calculate_inflections) {
    message("Calculating barcode inflection points by orig.ident...")
    mergedSeurat <- CalculateBarcodeInflections(mergedSeurat, group.column = 'orig.ident')
  }
  
  # Print summary statistics
  message("Merge summary:")
  message("  Total cells: ", ncol(mergedSeurat))
  message("  Total features: ", nrow(mergedSeurat))
  message("  Samples: ", paste(names(seurat_list), collapse = ", "))
  message("  Cells per sample:")
  print(table(mergedSeurat$sample_id))
  
  return(mergedSeurat)
}