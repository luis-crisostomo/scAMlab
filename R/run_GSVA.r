# Run GSVA on gene signatures
#' @export
run_GSVA <- function(seurat_obj, GeneSet, assay = "RNA", slot = "data", ...){
  
  # Load required libraries
  library(Seurat)
  library(GSVA)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  
  # Assuming you have:
  # - seurat_obj: your Seurat object with clusters in active.ident
  # - cell_type_signatures: a named list of gene signatures for each cell type
  
  # Step 1: Extract expression data from Seurat object
  # Get normalized counts (log-transformed)
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Step 2: Prepare gene signatures
  # Your cell_type_signatures should be a named list like:
  # cell_type_signatures <- list(
  #   "T_cells" = c("CD3D", "CD3E", "CD3G", "CD8A", "CD4"),
  #   "B_cells" = c("CD19", "MS4A1", "CD79A", "CD79B"),
  #   "NK_cells" = c("KLRD1", "NCAM1", "NKG7", "GNLY"),
  #   "Monocytes" = c("CD14", "FCGR3A", "CSF1R", "CD68"),
  #   "Dendritic_cells" = c("FCER1A", "CST3", "CLEC9A"),
  #   # ... add more cell types as needed
  # )
  
  # Step 3: Run GSVA
  # Convert Seurat expression matrix to regular matrix if needed
  expr_matrix <- as.matrix(expr_matrix)
  
  # Run GSVA with ssGSEA method (good for single-cell data)
  gsva_param <- ssgseaParam(
    expr_matrix,
    GeneSet,
    assay = NA_character_,
    annotation = NULL,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE,
    checkNA = c("auto", "yes", "no"),
    use ="na.rm"
  )
  
  gsva_scores <- gsva(gsva_param, 
                      verbose = TRUE,
                      ...)  # Adjust parallel processing as needed
  
  # Step 4: Add GSVA scores to Seurat metadata
  # Transpose scores matrix to have cells as rows
  gsva_scores_t <- t(gsva_scores)
  
  # Add to Seurat metadata
  for(cell_type in colnames(gsva_scores_t)) {
    seurat_obj <- AddMetaData(seurat_obj, 
                              metadata = gsva_scores_t[, cell_type], 
                              col.name = paste0("GSVA_", cell_type))
  }
  
  # Step 5: Calculate average scores per cluster
  cluster_scores <- data.frame()
  clusters <- levels(Idents(seurat_obj))
  
  for(cluster in clusters) {
    cluster_cells <- WhichCells(seurat_obj, idents = cluster)
    cluster_data <- gsva_scores_t[cluster_cells, , drop = FALSE]
    
    avg_scores <- colMeans(cluster_data)
    cluster_scores <- rbind(cluster_scores, avg_scores)
  }
  
  rownames(cluster_scores) <- clusters
  colnames(cluster_scores) <- colnames(gsva_scores_t)
  
  # Step 6: Assign cell types based on highest scores
  predicted_celltypes <- apply(cluster_scores, 1, function(x) {
    colnames(cluster_scores)[which.max(x)]
  })
  
  # Create a mapping dataframe
  cluster_annotation <- data.frame(
    cluster = names(predicted_celltypes),
    predicted_celltype = predicted_celltypes,
    max_score = apply(cluster_scores, 1, max)
  )
  
  print("Cluster annotations:")
  print(cluster_annotation)
  
  # Step 7: Visualizations
  
  # 7a. Heatmap of cluster scores
  pheatmap(cluster_scores, 
           scale = "column",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = "GSVA Scores per Cluster")
  
  # 7b. UMAP plots with GSVA scores
  # Example for plotting individual cell type scores
  # for(cell_type in names(cell_type_signatures)) {
  #   p <- FeaturePlot(seurat_obj, 
  #                    features = paste0("GSVA_", cell_type),
  #                    pt.size = 0.5) +
  #     ggtitle(paste("GSVA Score:", cell_type))
  #   print(p)
  # }
  
  # 7c. Add predicted cell types to Seurat object
  # Create a named vector for mapping clusters to cell types
  # celltype_mapping <- setNames(cluster_annotation$predicted_celltype, 
  #                              cluster_annotation$cluster)
  # 
  # # Add cell type annotations to metadata
  # seurat_obj$predicted_celltype <- celltype_mapping[as.character(Idents(seurat_obj))]
  # 
  # # Plot UMAP with predicted cell types
  # DimPlot(seurat_obj, 
  #         group.by = "predicted_celltype",
  #         label = TRUE,
  #         pt.size = 0.5) +
  #   ggtitle("Predicted Cell Types")
  # 
  # # Step 8: Quality control and refinement
  # 
  # # 8a. Check score distributions
  # # Boxplot of scores per cluster
  # library(tidyr)
  # scores_long <- cluster_scores %>%
  #   rownames_to_column("cluster") %>%
  #   pivot_longer(cols = -cluster, names_to = "celltype", values_to = "score")
  # 
  # ggplot(scores_long, aes(x = cluster, y = score, fill = celltype)) +
  #   geom_bar(stat = "identity", position = "dodge") +
  #   theme_minimal() +
  #   labs(title = "GSVA Scores by Cluster and Cell Type") +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # # 8b. Calculate confidence scores (difference between top 2 scores)
  # confidence_scores <- apply(cluster_scores, 1, function(x) {
  #   sorted_scores <- sort(x, decreasing = TRUE)
  #   sorted_scores[1] - sorted_scores[2]
  # })
  # 
  # cluster_annotation$confidence <- confidence_scores
  # 
  # # 8c. Flag low-confidence annotations
  # cluster_annotation$high_confidence <- cluster_annotation$confidence > 0.1  # Adjust threshold
  # 
  # print("Annotations with confidence scores:")
  # print(cluster_annotation)
  # 
  # # Step 9: Manual review and refinement
  # # You might want to manually review clusters with:
  # # - Low confidence scores
  # # - Unexpected cell type assignments
  # # - Check marker gene expression for validation
  # 
  # # Example: Check top marker genes for a specific cluster
  # cluster_to_check <- "0"  # Replace with cluster of interest
  # cluster_markers <- FindMarkers(seurat_obj, 
  #                                ident.1 = cluster_to_check,
  #                                min.pct = 0.25,
  #                                logfc.threshold = 0.25)
  # head(cluster_markers, 10)
  # 
  # # Step 10: Save results
  # # Save the annotated Seurat object
  # # saveRDS(seurat_obj, "annotated_seurat_object.rds")
  # 
  # # Save cluster annotations
  # # write.csv(cluster_annotation, "cluster_annotations.csv", row.names = FALSE)
  # 
  # print("Cell type annotation completed!")
  # print(paste("Total clusters:", length(clusters)))
  # print(paste("Unique predicted cell types:", length(unique(cluster_annotation$predicted_celltype))))
  # return(seurat_obj)
}
