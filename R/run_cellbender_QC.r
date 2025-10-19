#' Generate CellBender QC plots
#'
#' @description
#' Internal function that generates quality control plots comparing raw and
#' CellBender-corrected data for a single Seurat object. Creates six plots
#' showing cell probability metrics, feature/transcript comparisons, and
#' distribution analyses.
#'
#' @param seurat_obj A Seurat object containing both RAW (raw counts) and RNA
#'   (CellBender-corrected) assays. For Seurat v5, the RNA assay should have
#'   a counts layer.
#' @param sample Character string. Sample name used for file naming. Default: "Sample"
#' @param th.pct_diff_genes Numeric. Minimum percentage difference threshold
#'   for genes to include in Plot 2. Default: 25
#' @param th.feature Integer or NULL. Maximum number of genes to display in
#'   Plot 2. If NULL, shows all genes above th.pct_diff_genes. Default: NULL
#' @param th.pct_diff_cells Numeric. Minimum percentage difference threshold
#'   for cells to include in Plot 3. Default: 25
#' @param th.cell Integer or NULL. Maximum number of cells to display in
#'   Plot 3. If NULL, shows all cells above th.pct_diff_cells. Default: NULL
#' @param n_labels Integer. Number of top genes to label in Plot 2. Default: 10
#'
#' @return A patchwork object containing six ggplot2 plots arranged in a 3x2 grid
#'
#' @details
#' The function generates the following plots:
#' \itemize{
#'   \item Plot 1: Cell probability vs background fraction
#'   \item Plot 2: Raw vs CellBender features (genes) with top changed genes labeled
#'   \item Plot 3: Raw vs CellBender transcripts (cells)
#'   \item Plot 4: Cell probability vs droplet efficiency
#'   \item Plot 5: Genes detected per cell density comparison
#'   \item Plot 6: Cell UMI distribution density comparison
#' }
#'
#' The function also saves two CSV files:
#' \itemize{
#'   \item {sample}_Feature_stats.csv: Gene-level statistics
#'   \item {sample}_Cell_stats.csv: Cell-level statistics
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_abline geom_density
#'   scale_color_gradient scale_fill_manual scale_x_log10 scale_y_log10 labs
#'   theme_classic theme element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr %>%
#' @importFrom patchwork wrap_plots
#' @importFrom scales trans_format math_format
#'
#' @keywords internal
cellbender_QCplots <- function(seurat_obj, sample = "Sample",
                               th.pct_diff_genes = 25, th.feature = NULL,
                               th.pct_diff_cells = 25, th.cell = NULL,
                               n_labels = 10) {

  # Load required libraries
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
  library(scales)

  # Extract count matrices - Handle both Seurat v4 and v5
  # Check if RAW assay exists
  if (!"RAW" %in% names(seurat_obj@assays)) {
    stop("RAW assay not found in Seurat object")
  }

  # Get feature names from the Seurat object
  if (inherits(seurat_obj@assays$RAW, "Assay5")) {
    feature_names <- rownames(seurat_obj@assays$RAW)
  } else {
    feature_names <- rownames(seurat_obj@assays$RAW@counts)
  }

  # Get cell names
  if (inherits(seurat_obj@assays$RAW, "Assay5")) {
    cell_names <- colnames(seurat_obj@assays$RAW)
  } else {
    cell_names <- colnames(seurat_obj@assays$RAW@counts)
  }

  # Extract RAW counts
  if (inherits(seurat_obj@assays$RAW, "Assay5")) {
    # Seurat v5
    if ("counts" %in% names(seurat_obj@assays$RAW@layers)) {
      raw_counts <- seurat_obj@assays$RAW@layers$counts
    } else {
      stop("RAW assay does not have a 'counts' layer. Available layers: ",
           paste(names(seurat_obj@assays$RAW@layers), collapse = ", "))
    }
  } else {
    # Seurat v4
    raw_counts <- seurat_obj@assays$RAW@counts
  }

  # Ensure raw_counts is a matrix and set dimnames
  if (!is.matrix(raw_counts)) {
    raw_counts <- as.matrix(raw_counts)
  }
  rownames(raw_counts) <- feature_names
  colnames(raw_counts) <- cell_names

  # Extract CellBender counts from RNA assay
  if (!"RNA" %in% names(seurat_obj@assays)) {
    stop("RNA assay not found in Seurat object")
  }

  if (inherits(seurat_obj@assays$RNA, "Assay5")) {
    # Seurat v5
    if ("counts" %in% names(seurat_obj@assays$RNA@layers)) {
      cellbender_counts <- seurat_obj@assays$RNA@layers$counts
    } else {
      stop("RNA assay does not have a 'counts' layer. Available layers: ",
           paste(names(seurat_obj@assays$RNA@layers), collapse = ", "))
    }
  } else {
    # Seurat v4
    cellbender_counts <- seurat_obj@assays$RNA@counts
  }

  # Ensure cellbender_counts is a matrix and set dimnames
  if (!is.matrix(cellbender_counts)) {
    cellbender_counts <- as.matrix(cellbender_counts)
  }
  rownames(cellbender_counts) <- feature_names
  colnames(cellbender_counts) <- cell_names

  # Calculate metrics for each gene (feature)
  raw_features <- rowSums(raw_counts)
  cellbender_features <- rowSums(cellbender_counts)
  feature_diff <- raw_features - cellbender_features

  # Percentage difference with safe division
  pct_diff <- ifelse(raw_features == 0,
                     ifelse(cellbender_features == 0, 0, Inf),
                     (feature_diff / raw_features) * 100)

  # Create feature comparison table
  feature_comparison_table <- data.frame(
    Raw_features = raw_features,
    CellBender_features = cellbender_features,
    Feature_Diff = feature_diff,
    Pct_Diff = pct_diff,
    row.names = rownames(raw_counts)
  )

  feature_comparison_table <- feature_comparison_table[
    order(feature_comparison_table$Pct_Diff, decreasing = TRUE), ]

  # Save feature stats
  SaveData(feature_comparison_table,
           name = paste0(sample, "_Feature_stats"),
           type = "csv")

  print(head(feature_comparison_table, 10))

  # Calculate metrics for each cell
  raw_transcripts <- colSums(raw_counts)
  cellbender_transcripts <- colSums(cellbender_counts)
  transcript_diff <- raw_transcripts - cellbender_transcripts

  pct_diff <- ifelse(raw_transcripts == 0,
                     ifelse(cellbender_transcripts == 0, 0, Inf),
                     (transcript_diff / raw_transcripts) * 100)

  # Create cell comparison table
  cell_comparison_table <- data.frame(
    Raw_transcripts = raw_transcripts,
    CellBender_transcripts = cellbender_transcripts,
    Transcript_Diff = transcript_diff,
    Pct_Diff = pct_diff,
    row.names = colnames(raw_counts)
  )

  # Add cell metadata if available
  metadata_cols <- c("cell_probability", "background_fraction",
                     "cell_size", "droplet_efficiency")
  for (col in metadata_cols) {
    if (col %in% colnames(seurat_obj@meta.data)) {
      cell_comparison_table[[col]] <- seurat_obj@meta.data[[col]]
    }
  }

  cell_comparison_table <- cell_comparison_table[
    order(cell_comparison_table$Pct_Diff, decreasing = TRUE), ]

  # Save cell stats
  SaveData(cell_comparison_table,
           name = paste0(sample, "_Cell_stats"),
           type = "csv")

  print(head(cell_comparison_table, 10))

  # Plot 1: Cell probability vs background fraction
  plot1 <- ggplot(seurat_obj@meta.data,
                  aes(x = background_fraction, y = cell_probability)) +
    geom_point(alpha = 0.6, size = 0.8) +
    labs(title = "Cell Probability vs. Background Fraction",
         x = "Background Fraction",
         y = "Cell Probability") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  # Plot 2: Raw features vs CellBender features
  if (!is.null(th.feature)) {
    feature_plot_df <- head(
      subset(feature_comparison_table, Pct_Diff >= th.pct_diff_genes),
      th.feature)
  } else {
    feature_plot_df <- subset(feature_comparison_table,
                              Pct_Diff >= th.pct_diff_genes)
  }

  total_genes_remaining <- nrow(feature_comparison_table)
  genes_shown <- nrow(feature_plot_df)
  plot2_subtitle <- paste0("Showing ", genes_shown, " genes out of ",
                           total_genes_remaining, " remaining genes")

  feature_plot_df$gene_label <- ifelse(
    rownames(feature_plot_df) %in% rownames(head(feature_plot_df, n_labels)),
    rownames(feature_plot_df), "")

  plot2 <- ggplot(feature_plot_df,
                  aes(x = Raw_features, y = CellBender_features)) +
    geom_point(alpha = 0.6, size = 0.8, aes(color = Pct_Diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_gradient(low = "red", high = "blue",
                         name = "% transcripts\nremoved\n(Transcripts/Gene)") +
    scale_x_log10() +
    scale_y_log10() +
    geom_text_repel(aes(label = gene_label),
                    size = 3, max.overlaps = 15,
                    box.padding = 0.3, point.padding = 0.3,
                    segment.color = "grey50", segment.size = 0.2,
                    min.segment.length = 0) +
    labs(title = "Raw vs CellBender Features",
         subtitle = plot2_subtitle,
         x = "Raw Features", y = "CellBender Features") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.7) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  # Plot 3: Raw transcripts vs CellBender transcripts
  if (!is.null(th.cell)) {
    cell_plot_df <- head(
      subset(cell_comparison_table, Pct_Diff >= th.pct_diff_cells),
      th.cell)
  } else {
    cell_plot_df <- subset(cell_comparison_table,
                           Pct_Diff >= th.pct_diff_cells)
  }

  total_cells_remaining <- nrow(cell_comparison_table)
  cells_shown <- nrow(cell_plot_df)
  plot3_subtitle <- paste0("Showing ", cells_shown, " cells out of ",
                           total_cells_remaining, " remaining cells")

  plot3 <- ggplot(cell_plot_df,
                  aes(x = Raw_transcripts, y = CellBender_transcripts)) +
    geom_point(alpha = 0.6, size = 0.8, aes(color = Pct_Diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_gradient(low = "red", high = "blue",
                         name = "% transcripts\nremoved\n(Transcripts/Cell)") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Raw vs CellBender Transcripts",
         subtitle = plot3_subtitle,
         x = "Raw Transcripts", y = "CellBender Transcripts") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.7) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  # Plot 4: Cell probability vs Droplet efficiency
  plot4 <- ggplot(seurat_obj@meta.data,
                  aes(x = droplet_efficiency, y = cell_probability)) +
    geom_point(alpha = 0.6, size = 0.8) +
    labs(title = "Cell Probability vs. Droplet Efficiency",
         x = "Droplet Efficiency", y = "Cell Probability") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  # Plot 5: Number of genes detected per cell
  raw_genes_per_cell <- colSums(raw_counts > 0)
  cellbender_genes_per_cell <- colSums(cellbender_counts > 0)

  genes_per_cell_data <- data.frame(
    Genes = c(raw_genes_per_cell, cellbender_genes_per_cell),
    Data = rep(c("Raw", "CellBender"), each = length(raw_genes_per_cell))
  )

  plot5 <- ggplot(genes_per_cell_data, aes(x = Genes, fill = Data)) +
    geom_density(alpha = 0.7) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_fill_manual(values = c("Raw" = "lightcoral",
                                 "CellBender" = "lightblue")) +
    labs(title = "Genes Detected per Cell: Before vs After",
         x = "Number of Genes Detected", y = "Density") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  # Plot 6: Cell UMI distribution comparison
  cell_umi_data <- data.frame(
    UMIs = c(raw_transcripts, cellbender_transcripts),
    Data = rep(c("Raw", "CellBender"), each = length(raw_transcripts))
  )

  plot6 <- ggplot(cell_umi_data, aes(x = UMIs, fill = Data)) +
    geom_density(alpha = 0.7) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_fill_manual(values = c("Raw" = "lightcoral",
                                 "CellBender" = "lightblue")) +
    labs(title = "Cell UMI Distribution: Before vs After",
         x = "UMIs per Cell", y = "Density") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  # Combine plots
  combined_plot <- (plot1 | plot2 | plot3) / (plot4 | plot5 | plot6)
  print(combined_plot)

  SaveFigure(combined_plot, paste0(sample, "_cellbender_QC"),
             type = "tiff", width = 18, height = 10, res = 200)

  return(combined_plot)
}


#' Run CellBender QC Analysis
#'
#' @description
#' Generate quality control plots comparing raw and CellBender-corrected
#' single-cell RNA-seq data. Supports both single Seurat objects and lists
#' of Seurat objects with optional parallel processing.
#'
#' @param seurat_input A Seurat object or a list of Seurat objects. Each object
#'   must contain:
#'   \itemize{
#'     \item RAW assay: Raw count matrix
#'     \item RNA assay: CellBender-corrected count matrix
#'   }
#'   For Seurat v5 objects, the RNA assay must have a 'counts' layer.
#' @param sample_names Character vector. Names for each sample, used in file
#'   naming and plot titles. If NULL, defaults to "Sample" for single objects
#'   or "Sample_1", "Sample_2", etc. for lists. Default: NULL
#' @param param BiocParallelParam object for parallel processing (e.g.,
#'   MulticoreParam, SnowParam). If NULL, runs sequentially. Default: NULL
#' @param th.pct_diff_genes Numeric. Minimum percentage difference threshold
#'   for genes to include in the feature comparison plot. Default: 25
#' @param th.feature Integer or NULL. Maximum number of genes to display in
#'   the feature comparison plot. If NULL, shows all genes above
#'   th.pct_diff_genes threshold. Default: NULL
#' @param th.pct_diff_cells Numeric. Minimum percentage difference threshold
#'   for cells to include in the transcript comparison plot. Default: 25
#' @param th.cell Integer or NULL. Maximum number of cells to display in the
#'   transcript comparison plot. If NULL, shows all cells above
#'   th.pct_diff_cells threshold. Default: NULL
#' @param n_labels Integer. Number of top genes to label in the feature
#'   comparison plot. Default: 10
#'
#' @return
#' \itemize{
#'   \item For single Seurat object: Returns a patchwork plot object
#'   \item For list of Seurat objects: Returns a named list of patchwork plot objects
#' }
#'
#' @details
#' This function compares raw and CellBender-corrected data across multiple
#' dimensions:
#' \itemize{
#'   \item Gene-level changes: Total UMI counts per gene
#'   \item Cell-level changes: Total UMI counts per cell
#'   \item Detection rates: Number of genes detected per cell
#'   \item CellBender metrics: Cell probability, background fraction, droplet efficiency
#' }
#'
#' For each sample, the function saves:
#' \itemize{
#'   \item A TIFF image with 6 QC plots (18x10 inches, 200 DPI)
#'   \item A CSV file with gene-level statistics
#'   \item A CSV file with cell-level statistics
#' }
#'
#' @section Parallel Processing:
#' When processing multiple samples, you can enable parallel processing by
#' providing a BiocParallelParam object. For example:
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
#' # Single Seurat object
#' result <- run_cellbender_QC(
#'   seurat_input = seurat_obj,
#'   sample_names = "Sample1",
#'   th.pct_diff_genes = 30
#' )
#'
#' # Multiple Seurat objects with parallel processing
#' library(BiocParallel)
#' param <- MulticoreParam(workers = 4)
#'
#' results <- run_cellbender_QC(
#'   seurat_input = list(obj1, obj2, obj3),
#'   sample_names = c("Control", "Treatment1", "Treatment2"),
#'   param = param,
#'   th.pct_diff_genes = 25,
#'   n_labels = 15
#' )
#' }
#'
#' @seealso
#' \code{\link[BiocParallel]{BiocParallelParam}} for parallel processing options
#'
#' @export
run_cellbender_QC <- function(seurat_input,
                              sample_names = NULL,
                              param = NULL,
                              th.pct_diff_genes = 25,
                              th.feature = NULL,
                              th.pct_diff_cells = 25,
                              th.cell = NULL,
                              n_labels = 10) {

  # Check if input is a single Seurat object or a list
  if (inherits(seurat_input, "Seurat")) {
    # Single Seurat object
    cat("Processing single Seurat object...\n")
    sample_name <- ifelse(is.null(sample_names), "Sample", sample_names)
    result <- cellbender_QCplots(seurat_input,
                                 sample = sample_name,
                                 th.pct_diff_genes = th.pct_diff_genes,
                                 th.feature = th.feature,
                                 th.pct_diff_cells = th.pct_diff_cells,
                                 th.cell = th.cell,
                                 n_labels = n_labels)
    return(result)

  } else if (is.list(seurat_input)) {
    # Load BiocParallel for parallel processing
    library(BiocParallel)

    # List of Seurat objects
    n_objects <- length(seurat_input)
    cat(paste("Processing", n_objects, "Seurat objects...\n"))

    # Generate sample names if not provided
    if (is.null(sample_names)) {
      sample_names <- paste0("Sample_", seq_len(n_objects))
    } else if (length(sample_names) != n_objects) {
      warning("Length of sample_names doesn't match number of objects. Using default names.")
      sample_names <- paste0("Sample_", seq_len(n_objects))
    }

    if (!is.null(param)) {
      # Parallel processing
      cat("Running in parallel mode using BiocParallel...\n")

      worker_function <- function(i, seurat_list, sample_names,
                                  th.pct_diff_genes, th.feature,
                                  th.pct_diff_cells, th.cell, n_labels) {
        # Load required libraries in each worker
        library(ggplot2)
        library(ggrepel)
        library(dplyr)
        library(patchwork)
        library(scales)
        library(Seurat)

        # Load custom functions
        files.sources <- list.files("scripts/custom_functions",
                                    pattern = "\\.r$",
                                    full.names = TRUE)
        invisible(sapply(files.sources, source))

        cellbender_QCplots(seurat_list[[i]],
                           sample = sample_names[i],
                           th.pct_diff_genes = th.pct_diff_genes,
                           th.feature = th.feature,
                           th.pct_diff_cells = th.pct_diff_cells,
                           th.cell = th.cell,
                           n_labels = n_labels)
      }

      results <- bplapply(seq_len(n_objects),
                          worker_function,
                          seurat_list = seurat_input,
                          sample_names = sample_names,
                          th.pct_diff_genes = th.pct_diff_genes,
                          th.feature = th.feature,
                          th.pct_diff_cells = th.pct_diff_cells,
                          th.cell = th.cell,
                          n_labels = n_labels,
                          BPPARAM = param)

    } else {
      # Sequential processing
      cat("Running in sequential mode...\n")
      results <- list()
      for (i in seq_len(n_objects)) {
        cat(paste("Processing object", i, "of", n_objects, ":",
                  sample_names[i], "\n"))
        results[[i]] <- cellbender_QCplots(seurat_input[[i]],
                                           sample = sample_names[i],
                                           th.pct_diff_genes = th.pct_diff_genes,
                                           th.feature = th.feature,
                                           th.pct_diff_cells = th.pct_diff_cells,
                                           th.cell = th.cell,
                                           n_labels = n_labels)
      }
    }

    # Name the results list
    names(results) <- sample_names
    cat("Processing completed!\n")
    return(results)

  } else {
    stop("Input must be either a single Seurat object or a list of Seurat objects")
  }
}
