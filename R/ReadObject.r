#' Read an RDS object from the package’s data path
#'
#' @param name The base name (without .rds) of the file
#' @return The R object read from the file
#' @export
ReadObject <- function(name){
  full_path <- file.path(data_path, paste0(name, ".rds"))
  readRDS(full_path)
}
