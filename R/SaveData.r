#' Save a dataset to disk in multiple formats
#'
#' This function saves a data frame or other R object to a specified file path.
#' Supported formats are CSV, RDS, and XLSX (Excel).
#'
#' @param data The object to save. Typically a data frame, but can also be
#'   any R object if using \code{type = "rds"}.
#' @param name A character string with the file name (without extension).
#' @param type File type for saving. One of \code{"csv"}, \code{"rds"},
#'   or \code{"xlsx"}. Default is \code{"csv"}.
#'
#' @details
#' - The full save path is constructed from a global variable \code{data_path},
#'   combined with \code{name} and the chosen extension.
#' - If the directory does not exist, it is created automatically.
#' - For \code{type = "xlsx"}, the \pkg{openxlsx} package must be installed.
#'
#' @return Invisibly returns the full path to the saved file.
#'
#' @examples
#' \dontrun{
#' data_path <- "results"   # must be defined in your environment
#' SaveData(mtcars, "cars_data")         # saves as CSV
#' SaveData(mtcars, "cars_data", "rds")  # saves as RDS
#' SaveData(mtcars, "cars_data", "xlsx") # saves as Excel
#' }
#'
#' @export
SaveData <- function(data, name, type = "csv") {
  # Create full file path
  full_path <- file.path(data_path, paste0(name, ".", type))

  # Extract directory from the full path and create if it doesn't exist
  dir_path <- dirname(full_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  # Save the data
  if(type == "csv") {
    write.csv(data, full_path, row.names = FALSE)
  } else if(type == "rds") {
    saveRDS(data, full_path)
  } else if(type == "xlsx") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' needed for Excel files. Please install it.")
    }
    openxlsx::write.xlsx(data, full_path)
  } else {
    stop("Unsupported file type. Use 'csv', 'rds', or 'xlsx'.")
  }

  cat("Data saved to:", full_path, "\n")
}
