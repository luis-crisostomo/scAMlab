#' Low-Pass Gene and Cell Filtering for Single-Cell RNA-seq Data
#'
#' @description
#' Performs quality-based filtering on single-cell RNA-seq data by removing
#' low-expressed genes and optionally low-quality cells. Implements a two-step
#' filtering process: first removes genes based on expression criteria, then
#' optionally removes cells with insufficient transcripts after gene filtering.
#' Supports both single Seurat objects and lists of Seurat objects.
#'
#' @param seurat.obj A Seurat object or a list of Seurat objects to be filtered
#' @param assay Character string. Name of the assay to use for filtering.
#'   Default: "RNA"
#' @param min_cells Integer. Minimum number of cells in which a gene must be
#'   expressed (count > 0) to be retained. Genes expressed in fewer cells are
#'   removed. Default: 0 (no filtering based on cell count)
#' @param min_transcripts Integer. Minimum total transcript count across all
#'   cells for a gene to be retained. Genes with fewer total transcripts are
#'   removed. Default: 1 (removes genes with zero counts)
#' @param min_cell_transcripts Integer. Minimum total transcript count per cell
#'   after gene filtering. Cells with fewer transcripts are removed in the
#'   second filtering step. Set to 0 to skip cell filtering. Default: 0 (no
#'   cell filtering)
#' @param verbose Logical. If TRUE, prints detailed filtering statistics for
#'   each object. Default: TRUE
#'
#' @return
#' \itemize{
#'   \item For single Seurat object: Returns a filtered Seurat object
#'   \item For list of Seurat objects: Returns a named list of filtered Seurat
#'         objects with the same names as the input list
#' }
#'
#' @details
#' **Filtering Process:**
#'
#' The function implements a two-step filtering approach:
#'
#' \strong{Step 1: Gene Filtering}
#'
#' Genes are removed if they fail EITHER of these criteria:
#' \itemize{
#'   \item Expressed in fewer than \code{min_cells} cells (count > 0)
#'   \item Have fewer than \code{min_transcripts} total counts across all cells
#' }
#'
#' Both thresholds must be met for a gene to be retained. This dual-threshold
#' approach helps remove both:
#' \itemize{
#'   \item Sparsely expressed genes (expressed in very few cells)
#'   \item Very lowly expressed genes (low total counts even if in many cells)
#' }
#'
#' \strong{Step 2: Cell Filtering (optional)}
#'
#' After gene filtering, if \code{min_cell_transcripts > 0}, cells are removed
#' if they have fewer than \code{min_cell_transcripts} total counts in the
#' filtered gene set. This step is important because:
#' \itemize{
#'   \item Gene filtering may reduce some cells' transcript counts substantially
#'   \item Cells that fall below the threshold after gene filtering are likely
#'         low-quality or empty droplets
#' }
#'
#' @section Verbose Output:
#' When \code{verbose = TRUE}, the function prints detailed statistics:
#'
#' \strong{Gene Filtering Statistics:}
#' \itemize{
#'   \item Total genes in original object
#'   \item Genes expressed in < min_cells cells
#'   \item Genes with < min_transcripts total transcripts
#'   \item Breakdown of genes removed by each criterion
#'   \item Percentage of genes removed
#'   \item List of removed genes (up to 10)
#' }
#'
#' \strong{Cell Filtering Statistics (if applicable):}
#' \itemize{
#'   \item Total cells before filtering
#'   \item Cells with < min_cell_transcripts
#'   \item Percentage of cells removed
#'   \item Transcript count range in removed cells
#' }
#'
#' @section Threshold Guidelines:
#' Recommended starting thresholds depend on your data and experimental design:
#'
#' \strong{Conservative filtering (retain more features):}
#' \itemize{
#'   \item min_cells = 3-5
#'   \item min_transcripts = 5-10
#'   \item min_cell_transcripts = 200-500
#' }
#'
#' \strong{Moderate filtering (standard QC):}
#' \itemize{
#'   \item min_cells = 5-10
#'   \item min_transcripts = 10-20
#'   \item min_cell_transcripts = 500-1000
#' }
#'
#' \strong{Stringent filtering (high-quality cells only):}
#' \itemize{
#'   \item min_cells = 10-20
#'   \item min_transcripts = 20-50
#'   \item min_cell_transcripts = 1000-2000
#' }
#'
#' Note: Optimal thresholds vary by tissue type, sequencing depth, and
#' biological context. Examine QC plots before and after filtering.
#'
#' @section Safety Features:
#' The function includes safety checks to prevent data loss:
#' \itemize{
#'   \item If all genes would be removed, returns original object with a warning
#'   \item If all cells would be removed, skips cell filtering with a warning
#'   \item Validates that assay exists before filtering
#'   \item Checks for empty count matrices
#' }
#'
#' @section Processing Multiple Objects:
#' When processing a list of Seurat objects:
#' \itemize{
#'   \item Each object is filtered independently with the same thresholds
#'   \item Statistics are printed separately for each object (if verbose = TRUE)
#'   \item Original list names are preserved in the output
#'   \item If input list is unnamed, objects are named "Object_1", "Object_2", etc.
#' }
#'
#' @examples
#' \dontrun{
#' # Basic gene filtering: remove genes in < 3 cells
#' filtered_seurat <- lowpass_gene_filter(
#'   seurat.obj = seurat_obj,
#'   min_cells = 3
#' )
#'
#' # Filter genes AND cells
#' filtered_seurat <- lowpass_gene_filter(
#'   seurat.obj = seurat_obj,
#'   min_cells = 5,
#'   min_transcripts = 10,
#'   min_cell_transcripts = 500
#' )
#'
#' # Conservative filtering (retain more features)
#' filtered_seurat <- lowpass_gene_filter(
#'   seurat.obj = seurat_obj,
#'   min_cells = 3,
#'   min_transcripts = 5,
#'   min_cell_transcripts = 200
#' )
#'
#' # Stringent filtering (high-quality cells only)
#' filtered_seurat <- lowpass_gene_filter(
#'   seurat.obj = seurat_obj,
#'   min_cells = 20,
#'   min_transcripts = 50,
#'   min_cell_transcripts = 2000
#' )
#'
#' # Silent filtering (no verbose output)
#' filtered_seurat <- lowpass_gene_filter(
#'   seurat.obj = seurat_obj,
#'   min_cells = 10,
#'   min_transcripts = 20,
#'   verbose = FALSE
#' )
#'
#' # Process multiple samples
#' seurat_list <- list(
#'   Control = seurat_ctrl,
#'   Treatment1 = seurat_trt1,
#'   Treatment2 = seurat_trt2
#' )
#'
#' filtered_list <- lowpass_gene_filter(
#'   seurat.obj = seurat_list,
#'   min_cells = 5,
#'   min_transcripts = 10,
#'   min_cell_transcripts = 500
#' )
#'
#' # Access individual filtered objects
#' filtered_ctrl <- filtered_list$Control
#'
#' # Use different assay
#' filtered_seurat <- lowpass_gene_filter(
#'   seurat.obj = seurat_obj,
#'   assay = "SCT",
#'   min_cells = 5
#' )
#'
#' # Gene filtering only (no cell filtering)
#' filtered_seurat <- lowpass_gene_filter(
#'   seurat.obj = seurat_obj,
#'   min_cells = 10,
#'   min_transcripts = 20,
#'   min_cell_transcripts = 0  # Skip cell filtering
#' )
#' }
#'
#' @note
#' \itemize{
#'   \item This function modifies the Seurat object by subsetting features and
#'         cells. Other data (e.g., reductions, graphs) may need to be
#'         recalculated after filtering.
#'   \item The order of filtering (genes first, then cells) is intentional.
#'         Cell filtering based on post-gene-filtering transcript counts helps
#'         identify cells that were mainly expressing filtered genes.
#'   \item For very large datasets, consider filtering genes first with
#'         min_cells only, then examining cell quality separately.
#' }
#'
#' @seealso
#' \code{\link[Seurat]{subset}} for the underlying subsetting operation
#' \code{\link{run_babraham_qc}} for computing QC metrics before filtering
#' \code{\link{plot_QC_classic}} for visualizing QC metrics
#'
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom Matrix rowSums colSums
#'
#' @export
lowpass_gene_filter <- function(seurat.obj,
                                assay = "RNA",
                                min_cells = 0,    # Minimum number of cells expressing the gene
                                min_transcripts = 1, # Minimum total transcripts across all cells
                                min_cell_transcripts = 0, # Minimum transcripts per cell after gene filtering
                                verbose = TRUE) {

  # Helper function to process a single Seurat object
  filter_single_object <- function(obj, obj_name = NULL, assay, min_cells, min_transcripts, min_cell_transcripts, verbose) {
    # Validate input
    if (!inherits(obj, "Seurat")) {
      stop("Input must be a Seurat object")
    }

    # Check if assay exists
    if (!assay %in% names(obj@assays)) {
      stop(paste("Assay", assay, "not found in Seurat object"))
    }

    # Ensure we're working with the correct assay
    DefaultAssay(obj) <- assay

    # Get the counts matrix
    counts_matrix <- GetAssayData(obj, layer = "counts")

    # Validate counts matrix
    if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
      warning("Empty counts matrix detected")
      return(obj)
    }

    # Calculate row sums (sum of counts for each gene across all cells)
    gene_sums <- Matrix::rowSums(counts_matrix)

    # Calculate number of cells expressing each gene (cells with count > 0)
    cells_per_gene <- Matrix::rowSums(counts_matrix > 0)

    # Identify genes meeting both criteria:
    # 1. Expressed in at least min_cells number of cells
    # 2. Have at least min_transcripts total counts across all cells
    genes_to_keep <- (cells_per_gene >= min_cells) & (gene_sums >= min_transcripts)

    # Count genes that will be removed due to each criterion
    low_cells <- cells_per_gene < min_cells
    low_counts <- gene_sums < min_transcripts

    # Count how many genes will be removed
    n_total_genes <- nrow(counts_matrix)
    n_removed_cells <- sum(low_cells & !low_counts)  # Only due to cell threshold
    n_removed_counts <- sum(low_counts & !low_cells) # Only due to count threshold
    n_removed_both <- sum(low_cells & low_counts)    # Due to both thresholds
    n_removed_total <- sum(!genes_to_keep)

    percent_genes_removed <- round(n_removed_total / n_total_genes * 100, 2)

    if (verbose) {
      if (!is.null(obj_name)) {
        cat(paste0("\n=== PROCESSING: ", obj_name, " ===\n"))
      }
      cat("=== GENE FILTERING ===\n")
      cat(paste0("Total genes in original object: ", n_total_genes, "\n"))
      cat(paste0("Genes expressed in <", min_cells, " cells: ", sum(low_cells), "\n"))
      cat(paste0("Genes with <", min_transcripts, " total transcripts: ", sum(low_counts), "\n"))
      cat(paste0("   - Removed only due to cell threshold: ", n_removed_cells, "\n"))
      cat(paste0("   - Removed only due to transcript threshold: ", n_removed_counts, "\n"))
      cat(paste0("   - Removed due to both thresholds: ", n_removed_both, "\n"))
      cat(paste0("Total genes removed: ", n_removed_total, " (", percent_genes_removed, "%)\n"))
      cat(paste0("Genes remaining after filtering: ", sum(genes_to_keep), "\n"))

      # Optionally print some of the removed genes for verification
      if (n_removed_total > 0 && n_removed_total <= 10) {
        cat("Removed genes: ", paste(rownames(counts_matrix)[!genes_to_keep], collapse=", "), "\n")
      } else if (n_removed_total > 10) {
        cat("First 10 removed genes: ", paste(rownames(counts_matrix)[!genes_to_keep][1:10], collapse=", "), "...\n")
      }
    }

    # Check if all genes would be removed
    if (sum(genes_to_keep) == 0) {
      warning("All genes would be removed with current thresholds. Returning original object.")
      return(obj)
    }

    # Subset the Seurat object to keep only genes meeting criteria
    obj <- subset(obj, features = rownames(counts_matrix)[genes_to_keep])

    # Now filter cells if min_cell_transcripts > 0
    if (min_cell_transcripts > 0) {
      # Get the filtered counts matrix
      filtered_counts <- GetAssayData(obj, layer = "counts")

      # Calculate total transcripts per cell after gene filtering
      cell_transcripts <- Matrix::colSums(filtered_counts)

      # Identify cells to keep
      cells_to_keep <- cell_transcripts >= min_cell_transcripts

      # Count cells that will be removed
      n_total_cells <- ncol(filtered_counts)
      n_removed_cells_filter <- sum(!cells_to_keep)
      percent_cells_removed <- round(n_removed_cells_filter / n_total_cells * 100, 2)

      if (verbose) {
        cat("\n=== CELL FILTERING ===\n")
        cat(paste0("Total cells before cell filtering: ", n_total_cells, "\n"))
        cat(paste0("Cells with <", min_cell_transcripts, " transcripts: ", n_removed_cells_filter, "\n"))
        cat(paste0("Cells removed: ", n_removed_cells_filter, " (", percent_cells_removed, "%)\n"))
        cat(paste0("Cells remaining: ", sum(cells_to_keep), "\n"))

        # Show distribution of transcript counts in removed cells
        if (n_removed_cells_filter > 0) {
          removed_cell_counts <- cell_transcripts[!cells_to_keep]
          cat(paste0("Transcript count range in removed cells: ",
                     min(removed_cell_counts), " - ", max(removed_cell_counts), "\n"))
        }
      }

      # Check if all cells would be removed
      if (sum(cells_to_keep) == 0) {
        warning("All cells would be removed with current thresholds. Skipping cell filtering.")
      } else {
        # Subset cells
        obj <- subset(obj, cells = colnames(filtered_counts)[cells_to_keep])
      }
    }

    if (verbose) {
      cat(paste0("\n=== FINAL SUMMARY ===\n"))
      cat(paste0("Final object dimensions: ", nrow(obj), " genes x ", ncol(obj), " cells\n"))
    }

    return(obj)
  }

  # Main function logic: handle both single objects and lists
  if (inherits(seurat.obj, "Seurat")) {
    # Single Seurat object
    return(filter_single_object(seurat.obj, NULL, assay, min_cells, min_transcripts, min_cell_transcripts, verbose))

  } else if (is.list(seurat.obj)) {
    # List of Seurat objects
    if (length(seurat.obj) == 0) {
      stop("Empty list provided")
    }

    # Check if all elements are Seurat objects
    if (!all(sapply(seurat.obj, function(x) inherits(x, "Seurat")))) {
      stop("All elements in the list must be Seurat objects")
    }

    # Get names for the objects (use list names or create default names)
    obj_names <- names(seurat.obj)
    if (is.null(obj_names)) {
      obj_names <- paste0("Object_", seq_along(seurat.obj))
    }

    if (verbose) {
      cat(paste0("Processing ", length(seurat.obj), " Seurat objects...\n"))
    }

    # Process each object in the list
    filtered_objects <- list()
    for (i in seq_along(seurat.obj)) {
      filtered_objects[[i]] <- filter_single_object(
        seurat.obj[[i]],
        obj_names[i],
        assay,
        min_cells,
        min_transcripts,
        min_cell_transcripts,
        verbose
      )
    }

    # Preserve original names
    names(filtered_objects) <- obj_names

    if (verbose) {
      cat("\n=== ALL OBJECTS PROCESSED ===\n")
    }

    return(filtered_objects)

  } else {
    stop("Input must be either a Seurat object or a list of Seurat objects")
  }
}
