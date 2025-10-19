#' Data-Driven Variable Feature Selection for Single-Cell RNA-seq
#'
#' @description
#' Implements four different data-driven methods for selecting highly variable
#' features (HVFs) in single-cell RNA-seq data. Each method provides an
#' alternative approach to the standard Seurat variable feature selection,
#' allowing users to choose features based on variance explained, inflection
#' points, statistical thresholds, or PCA loadings. The function supports
#' multiple methods simultaneously for comparison.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data
#' @param method Integer vector specifying which method(s) to use. Options are:
#'   \itemize{
#'     \item 1: Elbow method on cumulative variance explained
#'     \item 2: Inflection point detection in variance distribution
#'     \item 3: Standardized Variance threshold on standardized variance
#'     \item 4: PCA-based feature selection with inflection point detection
#'   }
#'   Can be a single value (e.g., 3) or multiple values (e.g., c(1, 3)).
#'   Default: c(1, 3)
#' @param method_params Named list of parameters for each method. Each element
#'   should be named "method1", "method2", "method3", or "method4" and contain
#'   a list of method-specific parameters. See Details for available parameters.
#'   Default: list() (uses default parameters for all methods)
#' @param plot Logical. If TRUE, displays diagnostic plots for each method.
#'   Default: TRUE
#' @param save_plot Logical. If TRUE, saves the diagnostic plots as a TIFF file.
#'   Default: FALSE
#' @param plot_filename Character string. Base name for the saved plot file
#'   (will be converted to .tiff format). Default: "QC_varFeatures"
#' @param ... Additional arguments passed to \code{\link[Seurat]{NormalizeData}}
#'
#' @return A list containing results for each method used, with the following structure:
#' \itemize{
#'   \item \strong{methodX}: Results for each method, containing:
#'     \itemize{
#'       \item features: Character vector of selected feature names
#'       \item n_features: Number of features selected
#'       \item Method-specific metrics (see Details)
#'       \item parameters: Parameters used for this method
#'     }
#'   \item \strong{comparison}: (if multiple methods used) Comparison between methods:
#'     \itemize{
#'       \item common_features: Features selected by all methods
#'       \item method-specific differences
#'     }
#'   \item \strong{parameters}: Overall metadata about methods and parameters used
#'   \item \strong{plot_function}: Function to regenerate plots (if plot=TRUE)
#'   \item \strong{plot_dimensions}: Dimensions used for plots
#' }
#'
#' @details
#' **Method 1: Elbow Method on Cumulative Variance**
#'
#' Selects features based on cumulative variance explained, similar to choosing
#' principal components in PCA. Features are ranked by standardized variance,
#' and selection stops when cumulative variance reaches a threshold.
#'
#' Parameters (method_params$method1):
#' \itemize{
#'   \item max_features: Maximum number of features to consider (default: 8000)
#'   \item cumvar_threshold: Cumulative variance threshold (default: 0.8, i.e., 80%)
#' }
#'
#' **Method 2: Inflection Point in Variance Distribution**
#'
#' Identifies natural breakpoints in the variance distribution curve using
#' inflection point detection. This finds where the variance curve changes
#' curvature, indicating a transition between highly variable and less variable
#' genes. Multiple inflection points may be detected, requiring user selection
#' via interactive prompt.
#'
#' Parameters (method_params$method2):
#' \itemize{
#'   \item max_features: Maximum number of features to consider (default: 8000)
#' }
#'
#' Note: Requires the \code{inflection} package. User will be prompted to
#' select among detected inflection points.
#'
#' **Method 3: Standardized Variance Threshold**
#'
#' Selects features where standardized variance exceeds a fixed threshold.
#' This is the most straightforward approach, selecting all genes that show
#' variance above background levels.
#'
#' Parameters (method_params$method3):
#' \itemize{
#'   \item max_features: Maximum number of features to consider (default: 8000)
#'   \item variance_threshold: Minimum standardized variance (default: 1.3)
#' }
#'
#' **Method 4: PCA-Based Feature Selection**
#'
#' Uses PCA loadings to identify important features. First determines significant
#' PCs using an elbow plot, then calculates feature importance as the sum of
#' absolute loadings across these PCs. Inflection points in the feature
#' importance curve identify optimal cutoffs. Multiple inflection points may
#' be detected, requiring user selection via interactive prompt.
#'
#' Parameters (method_params$method4):
#' \itemize{
#'   \item max_features: Maximum number of features to consider (default: 8000)
#'   \item cumvar_threshold: Cumulative variance threshold for PC selection (default: 0.8)
#' }
#'
#' Note: Requires the \code{inflection} package. User will be prompted to
#' select among detected inflection points.
#'
#' @section Interactive Selection:
#' Methods 2 and 4 use inflection point detection, which may identify multiple
#' potential cutoffs. When this occurs:
#' \enumerate{
#'   \item Diagnostic plots are displayed showing all inflection points
#'   \item The function prints information about each inflection point
#'   \item User is prompted to select which inflection point to use
#'   \item Selection is stored in results$methodX$selected_features
#' }
#'
#' If only one inflection point is found, it is used automatically.
#'
#' @section Diagnostic Plots:
#' Each method generates diagnostic plots to help interpret results:
#' \itemize{
#'   \item \strong{Method 1}: Elbow plot showing cumulative variance with threshold line
#'   \item \strong{Method 2}: Variance distribution with inflection points marked
#'   \item \strong{Method 3}: Variance distribution colored by selection threshold
#'   \item \strong{Method 4}: Two plots - PCA elbow plot and feature importance
#'         with inflection points
#' }
#'
#' Plots are arranged in a grid layout (1-3 columns depending on number of methods).
#' When saved, plots are 5×5 inches per panel at 200 DPI in TIFF format.
#'
#' @section Method Comparison:
#' When multiple methods are used, the function automatically compares results:
#' \itemize{
#'   \item Prints number of features selected by each method
#'   \item Calculates intersection (common features across all methods)
#'   \item For 2 methods: identifies method-specific features
#'   \item Results stored in results$comparison
#' }
#'
#' @section Memory Efficiency:
#' The function creates a temporary DietSeurat object containing only the RNA
#' assay to reduce memory usage during computation. The original Seurat object
#' is not modified.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default methods (1 and 3)
#' hvf_results <- custom_VarFeatures(seurat_obj)
#'
#' # Access selected features from Method 1
#' selected_genes <- hvf_results$method1$features
#'
#' # Use a single method with custom parameters
#' hvf_results <- custom_VarFeatures(
#'   seurat_obj = seurat_obj,
#'   method = 3,
#'   method_params = list(
#'     method3 = list(
#'       max_features = 5000,
#'       variance_threshold = 1.5
#'     )
#'   )
#' )
#'
#' # Compare multiple methods
#' hvf_results <- custom_VarFeatures(
#'   seurat_obj = seurat_obj,
#'   method = c(1, 2, 3, 4),
#'   plot = TRUE,
#'   save_plot = TRUE,
#'   plot_filename = "HVF_comparison"
#' )
#'
#' # Check common features across methods
#' common_genes <- hvf_results$comparison$common_features
#' cat("Found", length(common_genes), "genes in common\n")
#'
#' # Use Method 1 with strict variance threshold
#' hvf_results <- custom_VarFeatures(
#'   seurat_obj = seurat_obj,
#'   method = 1,
#'   method_params = list(
#'     method1 = list(
#'       max_features = 10000,
#'       cumvar_threshold = 0.9  # 90% variance explained
#'     )
#'   )
#' )
#'
#' # Use Method 4 (PCA-based) with custom normalization
#' hvf_results <- custom_VarFeatures(
#'   seurat_obj = seurat_obj,
#'   method = 4,
#'   normalization.method = "LogNormalize",
#'   scale.factor = 10000
#' )
#'
#' # Apply selected features to original Seurat object
#' VariableFeatures(seurat_obj) <- hvf_results$method1$features
#' }
#'
#' @references
#' Inflection point detection: Christopoulos, D.T. (2016). "Developing methods
#' for identifying the inflection point of a convex/concave curve." arXiv:1206.5478
#'
#' @importFrom Seurat DietSeurat NormalizeData FindVariableFeatures
#'   VariableFeatures HVFInfo ScaleData RunPCA Loadings Idents
#' @importFrom inflection findiplist
#' @importFrom grDevices tiff dev.off rainbow
#' @importFrom graphics plot abline legend layout
#' @importFrom tools file_path_sans_ext
#' @importFrom stats na.omit
#'
#' @seealso
#' \code{\link[Seurat]{FindVariableFeatures}} for standard Seurat HVF selection
#' \code{\link[Seurat]{VariableFeatures}} for setting/getting variable features
#'
#' @export
custom_VarFeatures <- function(seurat_obj, method = c(1, 3),
                               method_params = list(),
                               plot = TRUE,
                               save_plot = FALSE,
                               plot_filename = "QC_varFeatures", ...) {

  # Validate inputs
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  valid_methods <- c(1, 2, 3, 4)
  if (!all(method %in% valid_methods)) {
    stop("Method must be 1, 2, 3, 4, or a vector containing these numbers")
  }

  # Set default parameters for each method
  default_params <- list(
    method1 = list(max_features = 8000, cumvar_threshold = 0.8),
    method2 = list(max_features = 8000),
    method3 = list(max_features = 8000, variance_threshold = 1.3),
    method4 = list(max_features = 8000, cumvar_threshold = 0.8)
  )

  # Merge user-provided parameters with defaults
  final_params <- default_params
  for (method_name in names(method_params)) {
    if (method_name %in% names(final_params)) {
      final_params[[method_name]] <- modifyList(final_params[[method_name]], method_params[[method_name]])
    } else {
      warning(paste("Unknown method parameter:", method_name, "- ignoring"))
    }
  }

  # Create a temporary copy using DietSeurat to reduce memory usage
  cat("Creating temporary Seurat object (RNA assay only) and normalizing data...\n")
  temp_obj <- DietSeurat(seurat_obj, assays = "RNA")

  # Normalize the data (pass additional arguments via ...)
  temp_obj <- NormalizeData(temp_obj, verbose = FALSE, ...)

  # Determine the maximum max_features needed across all methods
  max_features_needed <- max(sapply(final_params[paste0("method", method)], function(x) x$max_features %||% 8000))

  # Find variable features with the maximum number needed
  cat("Finding variable features with max_features =", max_features_needed, "...\n")
  temp_obj <- FindVariableFeatures(temp_obj, nfeatures = max_features_needed, verbose = FALSE)

  # Debug: Check if FindVariableFeatures worked
  var_features <- VariableFeatures(temp_obj)
  cat("Found", length(var_features), "variable features\n")

  if (length(var_features) == 0) {
    stop("FindVariableFeatures returned 0 features. Check your data and parameters.")
  }

  # Try to get HVF info and debug what's available
  hvf_info <- HVFInfo(temp_obj)
  cat("HVFInfo columns:", paste(colnames(hvf_info), collapse = ", "), "\n")
  cat("HVFInfo dimensions:", nrow(hvf_info), "x", ncol(hvf_info), "\n")

  # Check if the expected column exists with different possible names
  variance_col <- NULL
  possible_names <- c("vst.variance.standardized", "standardized_variance", "variance.standardized", "residual_variance")

  for (col_name in possible_names) {
    if (col_name %in% colnames(hvf_info)) {
      variance_col <- col_name
      cat("Using variance column:", variance_col, "\n")
      break
    }
  }

  if (is.null(variance_col)) {
    cat("Available columns in HVFInfo:\n")
    print(colnames(hvf_info))
    stop("Could not find variance column. Available columns printed above.")
  }

  results <- list()

  # Initialize plot storage
  plot_data <- list()

  # Method 1: Elbow method on variance explained
  if (1 %in% method) {
    cat("Running Method 1: Elbow method on cumulative variance...\n")

    # Get method-specific parameters
    params1 <- final_params$method1
    max_features <- params1$max_features
    cumvar_threshold <- params1$cumvar_threshold

    # Subset to variable features (limit to method-specific max_features)
    method1_features <- var_features[1:min(length(var_features), max_features)]
    var_data <- hvf_info[method1_features, , drop = FALSE]

    # Calculate cumulative variance explained
    variance_values <- var_data[[variance_col]]
    variance_values <- variance_values[!is.na(variance_values)]  # Remove any NAs

    if (length(variance_values) == 0) {
      stop("No valid variance values found. Check your data.")
    }

    cumvar <- cumsum(variance_values) / sum(variance_values)

    # Find elbow point (where cumulative variance reaches threshold)
    elbow_point <- which(cumvar >= cumvar_threshold)[1]
    if (is.na(elbow_point)) {
      elbow_point <- length(cumvar)
      warning(paste("Cumulative variance threshold of", cumvar_threshold, "not reached. Using all", length(cumvar), "features."))
    }

    results$method1 <- list(
      features = method1_features[1:elbow_point],
      n_features = elbow_point,
      cumvar_at_cutoff = cumvar[elbow_point],
      parameters = params1
    )

    cat("Method 1 selected", elbow_point, "features (", round(cumvar[elbow_point]*100, 1), "% variance explained)\n")

    # Store plot data
    plot_data$method1 <- list(
      x = 1:length(cumvar),
      y = cumvar,
      threshold = cumvar_threshold,
      cutoff = elbow_point,
      title = "Method 1: Elbow Plot",
      xlab = "Number of Features",
      ylab = "Cumulative Variance Explained"
    )
  }

  # Method 2: Inflection point in variance distribution
  if (2 %in% method) {
    cat("Running Method 2: Inflection point in variance distribution...\n")

    # Get method-specific parameters
    params2 <- final_params$method2
    max_features <- params2$max_features

    # Get all genes (not just variable features) but limit to max_features
    all_hvf_info <- HVFInfo(temp_obj)
    all_variance <- all_hvf_info[[variance_col]]
    names(all_variance) <- rownames(all_hvf_info)  # attach gene names
    all_variance <- all_variance[!is.na(all_variance)]

    # Sort by variance (descending) and limit to max_features
    sorted_variance <- sort(all_variance, decreasing = TRUE)
    sorted_variance <- sorted_variance[1:min(length(sorted_variance), max_features)]

    # Find inflection points (findiplist automatically tries both convex/concave)
    library(inflection)
    x <- 1:length(sorted_variance)
    y <- as.numeric(sorted_variance)

    # Single call to findiplist with performance improvement
    inflection_points <- findiplist(x, y, index = 0, doparallel = TRUE)
    inflection_points <- unique(na.omit(round(inflection_points)))

    # Create results for all inflection points
    inflection_results <- list()

    if (length(inflection_points) > 0) {
      for (i in seq_along(inflection_points)) {
        cutoff <- inflection_points[i]
        selected_features <- names(sorted_variance)[1:cutoff]

        inflection_results[[paste0("inflection_", i)]] <- list(
          features = selected_features,
          n_features = length(selected_features),
          cutoff_position = cutoff,
          variance_at_cutoff = sorted_variance[cutoff]
        )
      }
    } else {
      warning("No inflection points found for method 2.")
    }

    results$method2 <- list(
      inflection_points = inflection_points,
      all_inflections = inflection_results,
      parameters = params2
    )

    cat("Method 2 found", length(inflection_points), "inflection points at positions:", paste(inflection_points, collapse = ", "), "\n")

    # Store plot data
    plot_data$method2 <- list(
      x = 1:length(sorted_variance),
      y = sorted_variance,
      inflection_points = inflection_points,
      title = "Method 2: Inflection Point Detection",
      xlab = "Gene Rank",
      ylab = "Standardized Variance"
    )
  }

  # Method 3: Standardized Variance threshold
  if (3 %in% method) {
    cat("Running Method 3: Standardized Variance threshold...\n")

    # Get method-specific parameters
    params3 <- final_params$method3
    max_features <- params3$max_features
    variance_threshold <- params3$variance_threshold

    # Subset to variable features (limit to method-specific max_features)
    method3_features <- var_features[1:min(length(var_features), max_features)]
    var_data <- hvf_info[method3_features, , drop = FALSE]

    # Filter based on standardized variance > threshold
    variance_values <- var_data[[variance_col]]
    variance_values[is.na(variance_values)] <- 0  # Replace NAs with 0

    significant_mask <- variance_values > variance_threshold
    significant_features <- method3_features[significant_mask]

    results$method3 <- list(
      features = significant_features,
      n_features = length(significant_features),
      threshold_used = variance_threshold,
      parameters = params3
    )

    cat("Method 3 selected", length(significant_features), "features (standardized variance >", variance_threshold, ")\n")

    # Store plot data
    plot_data$method3 <- list(
      x = 1:length(method3_features),
      y = variance_values,
      significant_mask = significant_mask,
      threshold = variance_threshold,
      title = "Method 3: Variance Threshold",
      xlab = "Feature Rank",
      ylab = "Standardized Variance"
    )
  }

  # Method 4: PCA-based selection
  if (4 %in% method) {
    cat("Running Method 4: PCA-based feature selection with inflection point detection...\n")

    # Get method-specific parameters
    params4 <- final_params$method4
    max_features <- params4$max_features
    cumvar_threshold <- params4$cumvar_threshold

    # Use variable features for PCA (need to scale data first)
    # Limit to method-specific max_features
    method4_features <- var_features[1:min(length(var_features), max_features)]

    temp_obj <- ScaleData(temp_obj, features = method4_features, verbose = FALSE)
    temp_obj <- RunPCA(temp_obj, features = method4_features, npcs = 50, verbose = FALSE)

    # Determine significant PCs using elbow method
    pca_stdev <- temp_obj@reductions$pca@stdev
    pca_var_explained <- pca_stdev^2 / sum(pca_stdev^2)
    cumvar_pca <- cumsum(pca_var_explained)

    # Find elbow in PCA (where cumulative variance reaches 80% or flattens)
    pca_elbow <- which(cumvar_pca >= cumvar_threshold)[1]
    if (is.na(pca_elbow)) pca_elbow <- min(15, length(cumvar_pca))  # Default to 15 PCs max

    # Get loadings for significant PCs
    loadings <- Loadings(temp_obj, reduction = "pca")[, 1:pca_elbow, drop = FALSE]

    # Calculate feature importance as sum of absolute loadings across significant PCs
    feature_importance <- rowSums(abs(loadings))
    sorted_importance <- sort(feature_importance, decreasing = TRUE)

    # Find inflection points in the feature importance curve
    cat("Detecting inflection points in PCA loading importance curve...\n")
    library(inflection)

    x <- 1:length(sorted_importance)
    y <- as.numeric(sorted_importance)

    # Single call to findiplist with performance improvement
    inflection_points <- findiplist(x, y, index = 0, doparallel = TRUE)
    inflection_points <- unique(na.omit(round(inflection_points)))

    cat("Found", length(inflection_points), "inflection points at positions:", paste(inflection_points, collapse = ", "), "\n")

    # Create results for all inflection points
    inflection_results <- list()

    if (length(inflection_points) > 0) {
      for (i in seq_along(inflection_points)) {
        cutoff <- inflection_points[i]
        selected_features <- names(sorted_importance)[1:cutoff]

        inflection_results[[paste0("inflection_", i)]] <- list(
          features = selected_features,
          n_features = length(selected_features),
          cutoff_position = cutoff,
          importance_at_cutoff = sorted_importance[cutoff]
        )
      }
    } else {
      warning("No inflection points found for method 4.")
    }

    results$method4 <- list(
      inflection_points = inflection_points,
      all_inflections = inflection_results,
      pcs_used = pca_elbow,
      variance_explained_by_pcs = cumvar_pca[pca_elbow],
      feature_importance = sorted_importance,
      parameters = params4
    )

    cat("Method 4 completed using", pca_elbow, "PCs (",
        round(cumvar_pca[pca_elbow]*100, 1), "% variance explained)\n")

    # Store plot data for Method 4 (both PCA and feature importance plots)
    plot_data$method4 <- list(
      # PCA elbow plot data
      pca_x = 1:length(cumvar_pca),
      pca_y = cumvar_pca,
      pca_elbow = pca_elbow,
      # Feature importance plot data
      importance_x = 1:length(sorted_importance),
      importance_y = sorted_importance,
      inflection_points = inflection_points,
      title_pca = "Method 4: PCA Elbow Plot",
      title_importance = "Method 4: Feature Importance with Inflections"
    )
  }

  # Create and display plots FIRST (before user selection)
  if (plot || save_plot) {
    # Calculate layout based on number of methods
    n_methods <- length(method)
    n_plots <- n_methods + ifelse(4 %in% method, 1, 0)  # Method 4 has 2 plots

    # Determine layout
    if (n_plots == 1) {
      layout_matrix <- matrix(1, nrow = 1, ncol = 1)
    } else if (n_plots == 2) {
      layout_matrix <- matrix(1:2, nrow = 1, ncol = 2)
    } else if (n_plots <= 4) {
      layout_matrix <- matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)
    } else {
      # For more than 4 plots, use 3 columns
      nrows <- ceiling(n_plots / 3)
      layout_matrix <- matrix(1:(nrows * 3), nrow = nrows, ncol = 3, byrow = TRUE)
    }

    # Calculate plot dimensions based on layout (5x5 inches per plot)
    plot_width <- ncol(layout_matrix) * 5
    plot_height <- nrow(layout_matrix) * 5

    # Create the plotting function
    create_plots <- function() {
      layout(layout_matrix)

      # Plot Method 1
      if (1 %in% method && "method1" %in% names(plot_data)) {
        data <- plot_data$method1
        plot(data$x, data$y, type = "l",
             xlab = data$xlab, ylab = data$ylab, main = data$title)
        abline(h = data$threshold, col = "red", lty = 2)
        abline(v = data$cutoff, col = "blue", lty = 2)
        legend("bottomright",
               legend = c(paste("Threshold =", data$threshold),
                          paste("Selected =", data$cutoff, "features")),
               col = c("red", "blue"), lty = 2, cex = 0.7)
      }

      # Plot Method 2
      if (2 %in% method && "method2" %in% names(plot_data)) {
        data <- plot_data$method2
        plot(data$x, data$y, pch = 16, cex = 0.5,
             xlab = data$xlab, ylab = data$ylab, main = data$title)

        if (length(data$inflection_points) > 0) {
          colors <- rainbow(length(data$inflection_points))
          for (i in seq_along(data$inflection_points)) {
            abline(v = data$inflection_points[i], col = colors[i], lty = 2, lwd = 2)
          }
          legend("topright",
                 legend = paste("Inflection", 1:length(data$inflection_points), "=", data$inflection_points),
                 col = colors, lty = 2, cex = 0.6)
        }
      }

      # Plot Method 3
      if (3 %in% method && "method3" %in% names(plot_data)) {
        data <- plot_data$method3
        plot(data$x, data$y, pch = 16, cex = 0.5,
             xlab = data$xlab, ylab = data$ylab, main = data$title,
             col = ifelse(data$significant_mask, "red", "gray"))
        abline(h = data$threshold, col = "blue", lty = 2)
        legend("topright",
               legend = c(paste("Threshold =", data$threshold),
                          paste("Selected =", sum(data$significant_mask), "features")),
               col = "blue", lty = 2, cex = 0.7)
      }

      # Plot Method 4 (two plots)
      if (4 %in% method && "method4" %in% names(plot_data)) {
        data <- plot_data$method4

        # PCA elbow plot
        plot(data$pca_x, data$pca_y, type = "l",
             xlab = "Principal Component", ylab = "Cumulative Variance Explained",
             main = data$title_pca)
        abline(v = data$pca_elbow, col = "red", lty = 2)
        abline(h = data$pca_y[data$pca_elbow], col = "red", lty = 2)

        # Feature importance plot
        plot(data$importance_x, data$importance_y, pch = 16, cex = 0.5,
             xlab = "Feature Rank", ylab = "PCA Loading Importance",
             main = data$title_importance)

        if (length(data$inflection_points) > 0) {
          colors <- rainbow(length(data$inflection_points))
          for (i in seq_along(data$inflection_points)) {
            abline(v = data$inflection_points[i], col = colors[i], lty = 2, lwd = 2)
          }
          legend("topright",
                 legend = paste("Inflection", 1:length(data$inflection_points), "=", data$inflection_points),
                 col = colors, lty = 2, cex = 0.6)
        }
      }

      # Reset layout
      layout(1)
    }

    # Display plots if requested
    if (plot) {
      create_plots()
    }

    # Save plot if requested
    if (save_plot) {
      cat("Creating plot with dimensions:", plot_width, "x", plot_height, "inches\n")

      # Always save as TIFF for high quality
      if (!grepl("\\.(tiff|tif)$", plot_filename, ignore.case = TRUE)) {
        plot_filename <- paste0(tools::file_path_sans_ext(plot_filename), ".tiff")
        cat("Note: Filename changed to", plot_filename, "for TIFF format\n")
      }

      tiff(filename = paste0(fig_path,"/",plot_filename), width = plot_width, height = plot_height,
           units = "in", res = 200, compression = "lzw")
      create_plots()
      dev.off()

      cat("High-quality TIFF plot saved as:", plot_filename, "\n")
    }

    # Store the plot function in results for later use
    results$plot_function <- create_plots
    results$plot_dimensions <- c(width = plot_width, height = plot_height)
  }

  # Helper function to get features from a method
  get_method_features <- function(method_num, result_key) {
    if (!(result_key %in% names(results))) return(NULL)

    method_result <- results[[result_key]]

    # For methods that return multiple inflection points (2 and 4)
    if (method_num %in% c(2, 4) && "all_inflections" %in% names(method_result)) {
      if (length(method_result$all_inflections) == 0) {
        cat("Method", method_num, "found no inflection points.\n")
        return(NULL)
      } else if (length(method_result$all_inflections) == 1) {
        # Only one inflection point, use it automatically
        cat("Method", method_num, "has only one inflection point. Using it automatically.\n")
        return(method_result$all_inflections[[1]]$features)
      } else {
        # Multiple inflection points, ask user to select
        cat("\n=== INFLECTION POINT SELECTION ===\n")
        cat("Method", method_num, "found", length(method_result$all_inflections), "inflection points:\n")
        for (i in 1:length(method_result$all_inflections)) {
          infl_result <- method_result$all_inflections[[i]]
          cat("  ", i, ": Position", infl_result$cutoff_position, "(", infl_result$n_features, "features)\n")
        }
        cat("Please examine the plot(s) above to make an informed decision.\n")

        repeat {
          user_choice <- readline(paste("Select inflection point for Method", method_num, "(1-", length(method_result$all_inflections), "): "))
          choice_num <- as.numeric(user_choice)
          if (!is.na(choice_num) && choice_num >= 1 && choice_num <= length(method_result$all_inflections)) {
            selected_infl <- method_result$all_inflections[[choice_num]]
            cat("Using Method", method_num, "inflection point", choice_num, "(", selected_infl$n_features, "features)\n")
            return(selected_infl$features)
          } else {
            cat("Invalid choice. Please enter a number between 1 and", length(method_result$all_inflections), "\n")
          }
        }
      }
    } else {
      # Methods that return a single feature set (1 and 3)
      return(method_result$features)
    }
  }

  # NOW handle user selection for methods that require it (methods 2 and 4)
  # This happens AFTER plots are displayed
  for (m in method) {
    if (m %in% c(2, 4)) {
      result_key <- paste0("method", m)
      if (result_key %in% names(results)) {
        # Call get_method_features to trigger user prompt if needed
        selected_features <- get_method_features(m, result_key)

        # Store the selected features back in results for easy access
        if (!is.null(selected_features)) {
          results[[result_key]]$selected_features <- selected_features
          results[[result_key]]$selected_n_features <- length(selected_features)
        }
      }
    }
  }

  # If multiple methods requested, show comparison
  if (length(method) > 1) {
    cat("\nComparison of methods:\n")
    method_names <- c("1" = "Elbow", "2" = "Inflection", "3" = "Threshold", "4" = "PCA")

    for (m in method) {
      result_key <- paste0("method", m)
      if (result_key %in% names(results)) {
        # Handle methods that return multiple inflection points
        if (m %in% c(2, 4) && "selected_n_features" %in% names(results[[result_key]])) {
          cat("Method", m, "(", method_names[as.character(m)], "):", results[[result_key]]$selected_n_features, "features (selected)\n")
        } else if (m %in% c(2, 4) && "all_inflections" %in% names(results[[result_key]])) {
          n_inflections <- length(results[[result_key]]$all_inflections)
          cat("Method", m, "(", method_names[as.character(m)], "):", n_inflections, "inflection points found\n")
        } else {
          cat("Method", m, "(", method_names[as.character(m)], "):", results[[result_key]]$n_features, "features\n")
        }
      }
    }

    # Find intersections for all pairs of methods
    if (length(method) == 2) {
      m1 <- method[1]
      m2 <- method[2]
      m1_key <- paste0("method", m1)
      m2_key <- paste0("method", m2)

      # Get features (use selected_features if available, otherwise use get_method_features)
      features1 <- if ("selected_features" %in% names(results[[m1_key]])) {
        results[[m1_key]]$selected_features
      } else {
        get_method_features(m1, m1_key)
      }

      features2 <- if ("selected_features" %in% names(results[[m2_key]])) {
        results[[m2_key]]$selected_features
      } else {
        get_method_features(m2, m2_key)
      }

      if (!is.null(features1) && !is.null(features2)) {
        common_features <- intersect(features1, features2)
        cat("\nCommon features between methods:", length(common_features), "\n")

        results$comparison <- list(
          common_features = common_features,
          method1_only = setdiff(features1, features2),
          method2_only = setdiff(features2, features1)
        )
      }
    } else if (length(method) > 2) {
      # For more than 2 methods, find intersection of all
      all_features_list <- list()
      for (i in seq_along(method)) {
        m <- method[i]
        result_key <- paste0("method", m)

        # Get features (use selected_features if available, otherwise use get_method_features)
        features <- if ("selected_features" %in% names(results[[result_key]])) {
          results[[result_key]]$selected_features
        } else {
          get_method_features(m, result_key)
        }

        if (!is.null(features)) {
          all_features_list[[as.character(m)]] <- features
        }
      }

      if (length(all_features_list) > 1) {
        # Find intersection of all methods
        common_features <- all_features_list[[1]]
        for (i in 2:length(all_features_list)) {
          common_features <- intersect(common_features, all_features_list[[i]])
        }

        cat("\nCommon features across all methods:", length(common_features), "\n")
        results$comparison <- list(
          common_features = common_features,
          individual_method_features = all_features_list
        )
      }
    }
  }

  # Add metadata to results
  results$parameters <- list(
    methods_used = method,
    method_params = final_params[paste0("method", method)],
    variance_column_used = variance_col
  )

  return(results)
}
