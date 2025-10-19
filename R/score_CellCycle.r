#' Score Cell Cycle Phase for Single-Cell RNA-seq Data
#'
#' @description
#' Assigns cell cycle phase scores (S and G2M) to cells based on expression of
#' canonical cell cycle markers. For mouse data, automatically converts human
#' cell cycle genes to mouse orthologs using gprofiler2.
#'
#' @param seurat.obj A Seurat object
#' @param assay Character string. Name of the assay to use for scoring.
#'   Default: "RNA"
#' @param species Character string. Species of the data: "human" or "mouse".
#'   For mouse, human cell cycle genes are converted to mouse orthologs.
#'   Default: "mouse"
#'
#' @return Returns the Seurat object with added metadata columns:
#' \itemize{
#'   \item \strong{S.Score}: S phase score
#'   \item \strong{G2M.Score}: G2M phase score
#'   \item \strong{Phase}: Predicted cell cycle phase (G1, S, or G2M)
#'   \item \strong{CC.Difference}: Difference between S and G2M scores (S.Score - G2M.Score)
#' }
#'
#' @details
#' The function uses the updated 2019 cell cycle gene list from Seurat
#' (\code{cc.genes.updated.2019}). For mouse data, orthologous genes are
#' identified using gprofiler2's \code{gorth} function.
#'
#' Cell cycle phases are assigned based on relative scores:
#' \itemize{
#'   \item \strong{S phase}: High S.Score
#'   \item \strong{G2M phase}: High G2M.Score
#'   \item \strong{G1 phase}: Low scores for both S and G2M
#' }
#'
#' The \strong{CC.Difference} score is useful for:
#' \itemize{
#'   \item Regressing out cell cycle effects (negative values indicate G2M, positive indicate S)
#'   \item Identifying actively cycling vs. quiescent cells
#'   \item Visualizing cell cycle progression as a continuous variable
#' }
#'
#' If fewer than 5 genes are found for either S or G2M phase, the function
#' returns the object unchanged with a warning.
#'
#' @section Gene Availability:
#' The function reports how many cell cycle genes were found in the dataset.
#' Low gene recovery may indicate:
#' \itemize{
#'   \item Different gene naming conventions (e.g., Ensembl IDs instead of symbols)
#'   \item Aggressive gene filtering removed cell cycle markers
#'   \item Species mismatch
#' }
#'
#' @examples
#' \dontrun{
#' # Mouse data (default)
#' seurat_obj <- score_CellCycle(seurat_obj, species = "mouse")
#'
#' # Human data
#' seurat_obj <- score_CellCycle(seurat_obj, species = "human")
#'
#' # Visualize cell cycle phase
#' DimPlot(seurat_obj, group.by = "Phase")
#' FeaturePlot(seurat_obj, features = c("S.Score", "G2M.Score"))
#'
#' # Visualize as continuous variable
#' FeaturePlot(seurat_obj, features = "CC.Difference")
#'
#' # Regress out cell cycle effects during scaling
#' seurat_obj <- ScaleData(
#'   seurat_obj,
#'   vars.to.regress = c("S.Score", "G2M.Score")
#' )
#'
#' # Or regress out the difference score
#' seurat_obj <- ScaleData(
#'   seurat_obj,
#'   vars.to.regress = "CC.Difference"
#' )
#'
#' # Check distribution of cell cycle phases
#' table(seurat_obj$Phase)
#' }
#'
#' @references
#' Kolberg, L., Raudvere, U., Kuzmin, I. and  Vilo, J. (2020). gprofiler2 -- an R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler. F1000Research, 9(ELIXIR):709. DOI: https://doi.org/10.12688/f1000research.24956.2
#'
#' @importFrom Seurat CellCycleScoring
#' @importFrom gprofiler2 gorth
#'
#' @seealso
#' \code{\link[Seurat]{CellCycleScoring}} for the underlying scoring function
#'
#' @export
score_CellCycle <- function(seurat.obj, assay = 'RNA', species = "mouse"){

  # Helper function for cell cycle genes
  get_cell_cycle_genes <- function(species) {
    if (species == "human") {
      return(cc.genes.updated.2019)
    } else if (species == "mouse") {
      # Convert human genes to mouse gene symbols
      if (!requireNamespace("gprofiler2", quietly = TRUE)) {
        stop("gprofiler2 package is required for mouse cell cycle genes conversion")
      }

      s_genes <- gprofiler2::gorth(cc.genes.updated.2019$s.genes,
                                   source_organism = "hsapiens",
                                   target_organism = "mmusculus")$ortholog_name
      g2m_genes <- gprofiler2::gorth(cc.genes.updated.2019$g2m.genes,
                                     source_organism = "hsapiens",
                                     target_organism = "mmusculus")$ortholog_name

      return(list(s.genes = s_genes, g2m.genes = g2m_genes))
    }
  }

  # Get cell cycle genes based on species
  cc_genes <- get_cell_cycle_genes(species)
  s.genes <- cc_genes$s.genes
  g2m.genes <- cc_genes$g2m.genes

  # Check gene availability
  s.genes.in.data <- intersect(s.genes, rownames(seurat.obj))
  g2m.genes.in.data <- intersect(g2m.genes, rownames(seurat.obj))

  message(sprintf("S phase genes found: %d out of %d",
                  length(s.genes.in.data), length(s.genes)))
  message(sprintf("G2M phase genes found: %d out of %d",
                  length(g2m.genes.in.data), length(g2m.genes)))

  if (length(s.genes.in.data) < 5 || length(g2m.genes.in.data) < 5) {
    warning("Not enough cell cycle genes found. Skipping cell cycle scoring.")

    return(seurat.obj)

  } else {
    # Score cells for cell cycle phases
    message("Scoring cells for cell cycle phases...")
    seurat.obj <- CellCycleScoring(
      seurat.obj,
      s.features = s.genes.in.data,
      g2m.features = g2m.genes.in.data,
      assay = assay,
      set.ident = FALSE
    )

    # Calculate the difference between S and G2M scores
    message("Computing cell cycle difference score...")
    seurat.obj[["CC.Difference"]] <- seurat.obj$S.Score - seurat.obj$G2M.Score

    return(seurat.obj)
  }
}

