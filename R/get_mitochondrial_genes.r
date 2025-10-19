#' Retrieve mitochondrial genes using GO terms and Ensembl annotations
#'
#' This function queries the Ensembl BioMart database to retrieve genes
#' associated with mitochondrial compartments or processes, based on Gene
#' Ontology (GO) terms. Optionally, nuclear-encoded mitochondrial genes can
#' also be included.
#'
#' @param species Character. Ensembl species identifier, e.g. `"hsapiens"`
#'   for human, `"mmusculus"` for mouse (default: `"mmusculus"`).
#' @param include_nuclear Logical. Whether to include nuclear-encoded
#'   mitochondrial genes (default: `TRUE`).
#'
#' @details
#' The function relies on \pkg{biomaRt} to connect to the Ensembl database
#' and retrieve annotations. A predefined set of mitochondrial-related GO
#' terms is used, extended with energy metabolism terms if
#' `include_nuclear = TRUE`.
#'
#' @return A character vector of external gene names (symbols).
#'
#' @examples
#' \dontrun{
#' # Retrieve human mitochondrial genes
#' mito_genes <- get_mitochondrial_genes("hsapiens")
#'
#' # Retrieve mouse mitochondrial genes, including nuclear-encoded
#' mito_genes <- get_mitochondrial_genes("mmusculus", include_nuclear = TRUE)
#' }
#'
#' @export
get_mitochondrial_genes <- function(species = "mmusculus",
                                    include_nuclear = TRUE) {

  # Connect to Ensembl
  mart <- useMart("ensembl", dataset = paste0(species,"_gene_ensembl"))  # replace with species

  # Define GO terms for mitochondrial functions
  mito_go_terms <- c(
    "GO:0005739", # mitochondrion
    "GO:0005740", # mitochondrial envelope
    "GO:0005741", # mitochondrial outer membrane
    "GO:0005743", # mitochondrial inner membrane
    "GO:0005759", # mitochondrial matrix
    "GO:0005761", # mitochondrial ribosome
    "GO:0005762", # mitochondrial large ribosomal subunit
    "GO:0005763", # mitochondrial small ribosomal subunit
    "GO:0070469"  # respiratory chain
  )

  if (include_nuclear) {
    additional_terms <- c(
      "GO:0006120", # mitochondrial electron transport
      "GO:0042775", # mitochondrial ATP synthesis
      "GO:0006091", # generation of precursor metabolites and energy
      "GO:0055114"  # oxidation-reduction process (in mitochondrial context)
    )
    mito_go_terms <- c(mito_go_terms, additional_terms)
  }

  # Retrieve mitochondrial genes by GO
  mito_genes <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006"),
    filters = "go",
    values = mito_go_terms,
    mart = mart
  )

  unique(mito_genes$external_gene_name)

}
