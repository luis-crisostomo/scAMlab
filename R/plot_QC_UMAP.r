#' Generate Comprehensive QC UMAP Plots
#'
#' @description
#' Performs standard Seurat preprocessing (normalization, scaling, PCA,
#' clustering, UMAP) and generates a comprehensive set of quality control
#' visualizations arranged in a 4-row layout. This function is designed for
#' rapid QC assessment of single-cell RNA-seq data.
#'
#' @param seurat.obj A Seurat object containing single-cell RNA-seq data
#' @param nfeatures Either a numeric value specifying the number of variable
#'   features to identify (default: 8000), or a character vector of specific
#'   gene names to use as variable features
#' @param filename Character string. Base name for the output TIFF file
#'   (without extension). Default: "QC_UMAP"
#' @param genes_to_plot Character vector of length 2. Gene names to display
#'   in feature plots in the fourth row. Default: c("Ptprc", "Csf1r")
#' @param assay Character string. Name of the assay to use for analysis.
#'   Default: "RNA"
#'
#' @return Invisibly returns NULL. The function is called for its side effect
#'   of saving a multi-panel UMAP plot.
#'
#' @details
#' The function generates a 4-row visualization layout:
#'
#' **Row 1 (3 panels):**
#' \itemize{
#'   \item Cluster assignments with labels
#'   \item nCount_RNA distribution (discrete quantile bins: q5, q25, q50, q75, q95)
#'   \item nFeature_RNA distribution (discrete quantile bins: q5, q25, q50, q75, q95)
#' }
#'
#' **Row 2 (3 panels):**
#' \itemize{
#'   \item Mitochondrial percentage (percent.mt)
#'   \item Ribosomal percentage (percent.rb)
#'   \item Background fraction (if available) or contamination score
#' }
#'
#' **Row 3 (3 panels):**
#' \itemize{
#'   \item Novelty difference (log10GenesPerUMI difference from cluster median)
#'   \item Fast cluster assignments
#'   \item Doublet classification
#' }
#'
#' **Row 4 (3 panels):**
#' \itemize{
#'   \item Most common largest genes (genes with highest expression per cell,
#'         showing genes appearing in ≥500 cells)
#'   \item User-specified gene 1 expression
#'   \item User-specified gene 2 expression
#' }
#'
#' @section Processing Steps:
#' The function performs the following Seurat workflow:
#' \enumerate{
#'   \item Normalization (NormalizeData)
#'   \item Variable feature selection (FindVariableFeatures or user-defined)
#'   \item Scaling (ScaleData)
#'   \item PCA (RunPCA)
#'   \item Clustering (FindNeighbors + FindClusters at resolution 0.5)
#'   \item UMAP dimensionality reduction (RunUMAP, dims 1:20)
#' }
#'
#' @section Required Metadata:
#' The function expects the following metadata columns in the Seurat object:
#' \itemize{
#'   \item nCount_RNA, nFeature_RNA (automatically computed by Seurat)
#'   \item percent.mt, percent.rb (mitochondrial and ribosomal percentages)
#'   \item novelty_diff (difference from cluster median novelty score)
#'   \item fastcluster_clusters (preliminary cluster assignments)
#'   \item doublet_class (doublet detection results)
#'   \item largest_gene (gene with highest counts per cell)
#'   \item background_fraction (optional, from CellBender)
#'   \item contamination_score (optional, used if background_fraction is missing)
#' }
#'
#' @section Output:
#' Saves a TIFF file with dimensions 15×20 inches at 200 DPI resolution.
#' The file is saved using the SaveFigure function with the specified filename.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' plot_QC_UMAP(seurat_obj)
#'
#' # Custom number of variable features and output name
#' plot_QC_UMAP(
#'   seurat.obj = seurat_obj,
#'   nfeatures = 5000,
#'   filename = "Sample1_QC",
#'   genes_to_plot = c("Cd3e", "Ms4a1")
#' )
#'
#' # Using specific genes as variable features
#' custom_genes <- c("Gene1", "Gene2", "Gene3", ...)
#' plot_QC_UMAP(
#'   seurat.obj = seurat_obj,
#'   nfeatures = custom_genes,
#'   filename = "Custom_features_QC"
#' )
#' }
#'
#' @importFrom Seurat DefaultAssay NormalizeData FindVariableFeatures
#'   VariableFeatures ScaleData RunPCA FindNeighbors FindClusters RunUMAP
#'   DimPlot FeaturePlot Embeddings
#' @importFrom dplyr filter pull mutate case_when
#' @importFrom patchwork wrap_plots
#'
#' @seealso
#' \code{\link{run_cellbender_QC}} for CellBender-specific QC plots
#' \code{\link{plot_QC_classic}} for classical, Seurat-suggested QC plots
#' \code{\link{plot_QC_babraham}} for Babraham Institute QC plots (https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#QC)
#'
#' @export
plot_QC_UMAP <- function(seurat.obj, nfeatures = 8000, filename = "QC_UMAP", genes_to_plot = c("Ptprc", "Csf1r"), assay = "RNA"){

  # Set default assay
  DefaultAssay(seurat.obj) <- assay

  seurat.obj <- NormalizeData(seurat.obj)
  invisible(gc())

  if (length(nfeatures) == 1 & is.numeric(nfeatures)){

    seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = nfeatures)
    seurat.obj <- ScaleData(seurat.obj)
    invisible(gc())

  } else{
    VariableFeatures(seurat.obj) <- nfeatures
    seurat.obj <- ScaleData(seurat.obj)#, features = nfeatures)
    invisible(gc())
  }

  # Run PCA
  seurat.obj <- RunPCA(seurat.obj)
  invisible(gc())

  # clustering
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:20)
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)

  #UMAP plots
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:20)

  #Assemble plots
  UMAP.1.1 <- DimPlot(seurat.obj, reduction = "umap", repel = T, label = T)

  # Calculate quantiles for highlighting extreme values
  ncount_q5 <- quantile(seurat.obj$nCount_RNA, 0.05, na.rm = TRUE)
  ncount_q95 <- quantile(seurat.obj$nCount_RNA, 0.95, na.rm = TRUE)
  nfeature_q5 <- quantile(seurat.obj$nFeature_RNA, 0.05, na.rm = TRUE)
  nfeature_q95 <- quantile(seurat.obj$nFeature_RNA, 0.95, na.rm = TRUE)

  # Get UMAP coordinates and combine with metadata for extreme values
  umap_coords <- as.data.frame(Embeddings(seurat.obj, reduction = "umap"))
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

  # Combine UMAP coordinates with metadata
  plot_data <- cbind(seurat.obj@meta.data, umap_coords)

  # Calculate additional quantiles for binning
  ncount_q25 <- quantile(seurat.obj$nCount_RNA, 0.25, na.rm = TRUE)
  ncount_q50 <- quantile(seurat.obj$nCount_RNA, 0.50, na.rm = TRUE)
  ncount_q75 <- quantile(seurat.obj$nCount_RNA, 0.75, na.rm = TRUE)

  nfeature_q25 <- quantile(seurat.obj$nFeature_RNA, 0.25, na.rm = TRUE)
  nfeature_q50 <- quantile(seurat.obj$nFeature_RNA, 0.50, na.rm = TRUE)
  nfeature_q75 <- quantile(seurat.obj$nFeature_RNA, 0.75, na.rm = TRUE)

  # Define a discrete color palette (blue to red)
  qc_palette_count <- c(
    "blue",
    "lightblue",
    "lightgrey",
    "orange",
    "tomato",
    "darkred"
  )

  names(qc_palette_count) <- c(
    paste0("< ", round(ncount_q5, 0), " (q5)"),
    paste0("< ", round(ncount_q25, 0), " (q25)"),
    paste0("< ", round(ncount_q50, 0), " (q50)"),
    paste0("< ", round(ncount_q75, 0), " (q75)"),
    paste0("< ", round(ncount_q95, 0), " (q95)"),
    paste0("> ", round(ncount_q95, 0), " (q95)")
  )

  # Define a discrete color palette (blue to red)
  qc_palette_feature <- c(
    "blue",
    "lightblue",
    "lightgrey",
    "orange",
    "tomato",
    "darkred"
  )

  names(qc_palette_feature) <- c(
    paste0("< ", round(nfeature_q5, 0), " (q5)"),
    paste0("< ", round(nfeature_q25, 0), " (q25)"),
    paste0("< ", round(nfeature_q50, 0), " (q50)"),
    paste0("< ", round(nfeature_q75, 0), " (q75)"),
    paste0("< ", round(nfeature_q95, 0), " (q95)"),
    paste0("> ", round(nfeature_q95, 0), " (q95)")
  )

  # Add discrete bins to metadata
  seurat.obj$bin_nCount <- cut(
    seurat.obj$nCount_RNA,
    breaks = c(-Inf, ncount_q5, ncount_q25, ncount_q50, ncount_q75, ncount_q95, Inf),
    labels = names(qc_palette_count),  # use the labels directly!
    include.lowest = TRUE
  )

  seurat.obj$bin_nFeature <- cut(
    seurat.obj$nFeature_RNA,
    breaks = c(-Inf, nfeature_q5, nfeature_q25, nfeature_q50, nfeature_q75, nfeature_q95, Inf),
    labels = names(qc_palette_feature),  # use the labels directly!
    include.lowest = TRUE
  )

  # UMAP plots with discrete categories
  UMAP.1.2 <- DimPlot(seurat.obj, group.by = "bin_nCount", reduction = "umap") +
    scale_color_manual(values = qc_palette_count, name = "nCount_RNA") +
    ggtitle("nCounts_RNA")

  # UMAP.1.3: nFeature_RNA with quantile-based gradient
  UMAP.1.3 <- DimPlot(seurat.obj, group.by = "bin_nFeature", reduction = "umap") +
    scale_color_manual(values = qc_palette_feature, name = "nFeature_RNA") +
    ggtitle("nFeatures_RNA")

  # Issue 3: Replace cell_probability with background_fraction
  if(is.null(seurat.obj@meta.data[["background_fraction"]])){
    UMAP.2 <- FeaturePlot(seurat.obj, features = c("percent.mt", "percent.rb", "contamination_score"),
                          ncol = 3)
  } else {
    UMAP.2 <- FeaturePlot(seurat.obj, features = c("percent.mt", "percent.rb", "background_fraction"),
                          ncol = 3)
  }

  UMAP.3.1 <- (FeaturePlot(seurat.obj, features = "novelty_diff", reduction = "umap") +
                 scale_color_gradient2(low = "blue", mid = "lightgrey", high = "red",
                                       midpoint = 0,
                                       limits = c(-round(max(abs(seurat.obj$novelty_diff)), 1), round(max(abs(seurat.obj$novelty_diff)), 1)),  # Symmetric limits
                                       name = "novelty_diff"))
  UMAP.3.2 <- DimPlot(seurat.obj, reduction = "umap", repel = T, label = T, group.by = "fastcluster_clusters")
  UMAP.3.3 <- DimPlot(seurat.obj, reduction = "umap", repel = T, label = T, group.by = "doublet_class")

  UMAP.4.1 <-  ( {
    #plot the top most common "large" genes according to previous filter
    top_genes <- table(seurat.obj@meta.data$largest_gene) %>%
      as.data.frame() %>%
      filter(Freq >= 500) %>%
      pull(Var1) %>%
      as.character()

    # Add grouped column to metadata
    seurat.obj@meta.data <- seurat.obj@meta.data %>%
      mutate(
        largest_gene_plot = case_when(
          largest_gene %in% top_genes ~ as.character(largest_gene),
          TRUE ~ "Other"
        ),
        largest_gene_plot = factor(largest_gene_plot,
                                   levels = c(top_genes, "Other"))
      )

    DimPlot(seurat.obj, reduction = "umap", repel = T, label = T, group.by = "largest_gene_plot")
  }
  # FeaturePlot(seurat.obj, features = names(head(sort(table(seurat.obj@meta.data$largest_gene), decreasing = T), 15)),
  #             reduction = "umap")
  )

  # Issue 3: Use user-defined genes for plotting
  UMAP.4.2 <- FeaturePlot(seurat.obj, features = genes_to_plot[1], reduction = "umap")
  UMAP.4.3 <- FeaturePlot(seurat.obj, features = genes_to_plot[2], reduction = "umap")

  # Issue 1: Use wrap_plots() to ensure proper 3-column layout for third row
  library(patchwork)

  row1 <- UMAP.1.1 + UMAP.1.2 + UMAP.1.3
  row2 <- UMAP.2
  row3 <- wrap_plots(UMAP.3.1, UMAP.3.2, UMAP.3.3, ncol = 3)
  row4 <- UMAP.4.1 + UMAP.4.2 + UMAP.4.3

  combine.plot <- row1 / row2 / row3 / row4
  #combine.plot

  SaveFigure(combine.plot, name = filename, type = "tiff", width = 15, height = 20, res = 200)

}
