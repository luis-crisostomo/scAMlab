#' @export
SaveFigure <- function(plots, name, type = "tiff", width, height, res = 200) {
  # Create full file path
  full_path <- file.path(fig_path, paste0(name, ".", type))
  
  # Extract directory from the full path and create if it doesn't exist
  dir_path <- dirname(full_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save the plot
  if(type == "tiff") {
    tiff(full_path, width = width, height = height, units = "in", res = res)
  } else if(type == "pdf") {
    pdf(full_path, width = width, height = height)
  } else if(type == "png") {
    png(full_path, width = width, height = height, units = "in", res = res)
  } else {
    stop("Unsupported file type. Use 'tiff', 'pdf', or 'png'.")
  }
  
  print(plots)
  dev.off()
  
  cat("Figure saved to:", full_path, "\n")
}