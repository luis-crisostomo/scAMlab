#' Retrieve ribosomal genes using GO terms and Ensembl annotations
#'
#' This function queries the Ensembl BioMart database to retrieve genes
#' associated with ribosomal structure or translation, based on Gene
#' Ontology (GO) terms. Optionally, more specific ribosomal subunit-related
#' terms can be included.
#'
#' @param species Character. Ensembl species identifier, e.g. `"hsapiens"`
#'   for human, `"mmusculus"` for mouse (default: `"mmusculus"`).
#' @param include_subunits Logical. Whether to include specific ribosomal
#'   subunit and ribonucleoprotein-related terms (default: `TRUE`).
#'
#' @details
#' The function relies on \pkg{biomaRt} to connect to the Ensembl database
#' and retrieve annotations. A predefined set of ribosome-related GO terms
#' is used, extended with cytosolic and ribonucleoprotein terms if
#' `include_subunits = TRUE`.
#'
#' @return A character vector of external gene names (symbols).
#'
#' @examples
#' \dontrun{
#' # Retrieve human ribosomal genes
#' ribo_genes <- get_ribosomal_genes("hsapiens")
#'
#' # Retrieve mouse ribosomal genes with additional subunit terms
#' ribo_genes <- get_ribosomal_genes("mmusculus", include_subunits = TRUE)
#' }
#'
#' @export
get_ribosomal_genes <- function(species = "mmusculus", include_subunits = TRUE){

  # Connect to Ensembl
  mart <- useMart("ensembl", dataset = paste0(species,"_gene_ensembl"))  # replace with species

  # Define GO terms for ribosome-related functions
  ribosome_go_terms <- c(
    "GO:0005840", # ribosome
    "GO:0003735", # structural constituent of ribosome
    "GO:0006412", # translation
    "GO:0022625", # cytosolic large ribosomal subunit
    "GO:0022627", # cytosolic small ribosomal subunit
    "GO:0000315", # organellar large ribosomal subunit
    "GO:0000314"  # organellar small ribosomal subunit
  )

  if (include_subunits) {
    additional_terms <- c(
      "GO:0005829", # cytosol (for cytosolic ribosomes)
      "GO:0030529", # ribonucleoprotein complex
      "GO:0044391"  # ribosomal subunit
    )
    ribosome_go_terms <- c(ribosome_go_terms, additional_terms)
  }

  # Retrieve ribosomal genes by GO
  ribosomal_genes <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006"),
    filters = "go",
    values = ribosome_go_terms,
    mart = mart
  )

  unique(ribosomal_genes$external_gene_name)

}
