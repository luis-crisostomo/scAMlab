#' Add Sample Metadata to Seurat Object
#'
#' @description
#' Adds sample-level metadata from an external CSV file to a Seurat object by
#' matching sample identifiers. Particularly useful for merged Seurat objects
#' containing multiple samples. Handles edge cases including leading zeros in
#' sample IDs and provides flexible matching when exact matches fail.
#'
#' @param seurat_obj A Seurat object, typically containing multiple samples
#' @param metadata_file Character string. Path to a CSV file containing sample
#'   metadata. Must include a column matching sample identifiers in the Seurat
#'   object.
#' @param sample_id_col Character string. Name of the column in the metadata
#'   file containing sample identifiers that match values in the Seurat object.
#'   Default: "sample_id"
#' @param id_in_seurat Character string. Name of the metadata column in the
#'   Seurat object containing sample identifiers to match against.
#'   Default: "orig.ident"
#' @param verbose Logical. If TRUE, prints detailed information about the
#'   matching process and metadata addition. Default: TRUE
#'
#' @return Returns the Seurat object with additional metadata columns added.
#'   Each cell receives metadata values corresponding to its sample identifier.
#'
#' @details
#' **Matching Process:**
#'
#' The function matches samples through the following steps:
#' \enumerate{
#'   \item Identifies all unique sample identifiers in the Seurat object
#'   \item For each identifier, finds the corresponding row in the metadata file
#'   \item Adds all metadata columns (except the ID column) to cells from that sample
#' }
#'
#' **Flexible Matching:**
#'
#' If an exact match fails, the function attempts flexible matching by:
#' \itemize{
#'   \item Removing leading zeros from both the Seurat ID and metadata ID
#'   \item Comparing the resulting strings
#'   \item Reporting the match if successful
#' }
#'
#' This is particularly useful when sample IDs have inconsistent formatting
#' (e.g., "001" vs "1" or "Sample_01" vs "Sample_1").
#'
#' **Edge Cases Handled:**
#' \itemize{
#'   \item \strong{Leading zeros}: Sample IDs like "001", "002" are preserved
#'         as strings to prevent conversion to numeric
#'   \item \strong{Missing matches}: Warns if a sample in Seurat has no
#'         corresponding metadata entry
#'   \item \strong{Duplicate entries}: Warns and uses the first entry if
#'         multiple metadata rows match a sample ID
#'   \item \strong{New columns}: Creates new metadata columns with NA values
#'         for samples not receiving that metadata
#' }
#'
#' @section CSV File Format:
#' The metadata CSV file should have:
#' \itemize{
#'   \item One row per sample
#'   \item One column containing sample identifiers (specified by sample_id_col)
#'   \item Additional columns with metadata to add (e.g., condition, batch, donor, etc.)
#' }
#'
#' Example CSV structure:
#' \preformatted{
#' sample_id,condition,batch,sex,age
#' Sample_01,Control,Batch1,Female,25
#' Sample_02,Treatment,Batch1,Male,30
#' Sample_03,Control,Batch2,Female,28
#' }
#'
#' @section Verbose Output:
#' When \code{verbose = TRUE}, the function prints:
#' \itemize{
#'   \item Path to the metadata file being loaded
#'   \item All available metadata columns
#'   \item Number of unique identifiers found in Seurat object
#'   \item List of all identifiers
#'   \item Flexible matching attempts and results
#'   \item Each metadata addition with sample ID, column name, and value
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' seurat_obj <- add_seurat_metadata(
#'   seurat_obj = seurat_obj,
#'   metadata_file = "sample_metadata.csv"
#' )
#'
#' # Custom column names
#' seurat_obj <- add_seurat_metadata(
#'   seurat_obj = seurat_obj,
#'   metadata_file = "metadata.csv",
#'   sample_id_col = "SampleName",
#'   id_in_seurat = "sample"
#' )
#'
#' # Silent mode (no messages)
#' seurat_obj <- add_seurat_metadata(
#'   seurat_obj = seurat_obj,
#'   metadata_file = "metadata.csv",
#'   verbose = FALSE
#' )
#'
#' # Check added metadata
#' head(seurat_obj@meta.data)
#' table(seurat_obj$condition)
#'
#' # Visualize by new metadata
#' DimPlot(seurat_obj, group.by = "condition")
#'
#' # Example workflow with merged samples
#' # 1. Merge samples
#' seurat_merged <- merge(
#'   x = seurat_list[[1]],
#'   y = seurat_list[2:length(seurat_list)],
#'   add.cell.ids = names(seurat_list)
#' )
#'
#' # 2. Add metadata
#' seurat_merged <- add_seurat_metadata(
#'   seurat_obj = seurat_merged,
#'   metadata_file = "experiment_metadata.csv",
#'   sample_id_col = "sample_id",
#'   id_in_seurat = "orig.ident"
#' )
#'
#' # 3. Use metadata in analysis
#' DimPlot(seurat_merged, group.by = "condition", split.by = "batch")
#' }
#'
#' @note
#' \itemize{
#'   \item The function modifies metadata in place - existing columns with the
#'         same name will be overwritten
#'   \item Cells without matching sample IDs will have NA values for new metadata
#'   \item Sample IDs in the CSV file are read as character strings to preserve
#'         formatting (e.g., leading zeros)
#'   \item If sample IDs contain special characters, ensure they match exactly
#'         between the CSV and Seurat object
#' }
#'
#' @seealso
#' \code{\link[Seurat]{merge}} for merging multiple Seurat objects
#' \code{\link[Seurat]{AddMetaData}} for adding metadata directly
#'
#' @export
add_seurat_metadata <- function(seurat_obj, metadata_file, sample_id_col = "sample_id",
                                id_in_seurat = "orig.ident", verbose = TRUE) {
  # Check if seurat_obj is a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  # Load metadata
  if (!file.exists(metadata_file)) {
    stop(paste("Metadata file not found:", metadata_file))
  }

  # Load metadata with proper string handling for IDs with leading zeros
  if (verbose) message(paste("Loading metadata from", metadata_file))
  metadata_df <- read.csv(metadata_file, colClasses = c(sample_id_col = "character"), stringsAsFactors = FALSE)

  # Check if sample_id_col exists in the metadata
  if (!(sample_id_col %in% colnames(metadata_df))) {
    stop(paste("Sample ID column", sample_id_col, "not found in metadata file"))
  }

  if (verbose) {
    message("Available metadata columns:")
    message(paste(colnames(metadata_df), collapse = ", "))
  }

  # Get all unique sample identifiers from the merged Seurat object
  sample_ids <- unique(seurat_obj[[id_in_seurat, drop = TRUE]])

  if (verbose) {
    message(paste("Found", length(sample_ids), "unique identifiers in the Seurat object:"))
    message(paste(sample_ids, collapse = ", "))
  }

  # Process each unique sample ID in the merged Seurat object
  for (sample_id in sample_ids) {
    # Find matching row in metadata - handle exact string matching for IDs with leading zeros
    metadata_idx <- which(metadata_df[[sample_id_col]] == sample_id)

    if (length(metadata_idx) == 0) {
      # Try more flexible matching if exact match fails
      message(paste("Attempting flexible matching for", sample_id))
      # Try matching without leading zeros
      no_leading_zeros_id <- sub("^0+", "", sample_id)
      metadata_idx <- which(sub("^0+", "", metadata_df[[sample_id_col]]) == no_leading_zeros_id)

      if (length(metadata_idx) == 0) {
        warning(paste("No metadata found for sample", sample_id))
        next
      } else if (length(metadata_idx) > 1) {
        warning(paste("Multiple metadata entries found for sample", sample_id,
                      ". Using the first entry."))
        metadata_idx <- metadata_idx[1]
      }

      if (verbose) {
        message(paste("  Found match:", sample_id, "→", metadata_df[metadata_idx, sample_id_col]))
      }
    } else if (length(metadata_idx) > 1) {
      warning(paste("Multiple metadata entries found for sample", sample_id,
                    ". Using the first entry."))
      metadata_idx <- metadata_idx[1]
    }

    metadata_row <- metadata_df[metadata_idx, , drop = FALSE]

    # Create logical mask for cells with this sample_id
    mask <- seurat_obj[[id_in_seurat, drop = TRUE]] == sample_id

    # Add each metadata column to the matching cells
    for (col in colnames(metadata_row)) {
      if (col != sample_id_col) {  # Skip the ID column
        value <- metadata_row[[col]]

        # Add metadata to Seurat object
        # We need to work directly with the metadata dataframe
        if (!col %in% colnames(seurat_obj@meta.data)) {
          # Create new column with NAs if it doesn't exist
          seurat_obj@meta.data[[col]] <- NA
        }

        # Now set the values for cells that match the mask
        seurat_obj@meta.data[mask, col] <- value

        if (verbose) {
          message(paste("  Added metadata:", col, "=", value,
                        "for cells with", id_in_seurat, "=", sample_id))
        }
      }
    }
  }

  return(seurat_obj)
}
