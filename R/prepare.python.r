# ---------- Helper function to find conda ----------
find_conda_path <- function() {
  sys_name <- Sys.info()[["sysname"]]
  is_windows <- sys_name == "Windows"

  if (is_windows) {
    conda_home <- Sys.getenv("CONDA_HOME")
    if (conda_home == "") stop("CONDA_HOME environment variable is not set.")
    conda_home <- normalizePath(conda_home, "/")

    # Find mamba executable (can be in Scripts or Library/bin)
    mamba_paths <- c(
      file.path(conda_home, "Scripts", "mamba.exe"),
      file.path(conda_home, "Library", "bin", "mamba.exe"),
      file.path(conda_home, "Scripts", "mamba"),
      file.path(conda_home, "Library", "bin", "mamba")
    )

    mamba_exe <- NULL
    for (path in mamba_paths) {
      if (file.exists(path)) {
        mamba_exe <- path
        break
      }
    }

    return(list(
      conda_base = conda_home,
      conda_exe = file.path(conda_home, "Scripts", "conda.exe"),
      mamba_exe = mamba_exe
    ))
  }

  # For Linux/macOS, try multiple approaches to find conda
  conda_paths <- c(
    "~/miniconda3",
    "~/anaconda3",
    "/opt/conda",
    "/opt/miniconda3",
    "/opt/anaconda3",
    "/usr/local/miniconda3",
    "/usr/local/anaconda3"
  )

  # First, try common installation paths
  for (path in conda_paths) {
    expanded_path <- path.expand(path)
    conda_exe <- file.path(expanded_path, "bin", "conda")
    if (file.exists(conda_exe)) {
      return(list(
        conda_base = expanded_path,
        conda_exe = conda_exe,
        mamba_exe = file.path(expanded_path, "bin", "mamba")
      ))
    }
  }

  # Try to find conda in user's shell PATH by checking common shell rc files
  shell_files <- c("~/.bashrc", "~/.bash_profile", "~/.zshrc", "~/.profile")
  for (shell_file in shell_files) {
    expanded_shell_file <- path.expand(shell_file)
    if (file.exists(expanded_shell_file)) {
      content <- readLines(expanded_shell_file, warn = FALSE)

      for (line in content) {
        if (grepl("^\\s*#", line) || grepl("^\\s*$", line)) next

        # Pattern: source .../etc/profile.d/conda.sh
        if (grepl("source.*conda\\.sh", line) || grepl("\\. .*conda\\.sh", line)) {
          path_match <- regmatches(line, regexpr("['\"]?([^'\"\\s]+)/etc/profile\\.d/conda\\.sh['\"]?", line, perl = TRUE))
          if (length(path_match) > 0) {
            conda_base <- gsub("['\"]", "", path_match[1])
            conda_base <- gsub("/etc/profile\\.d/conda\\.sh.*$", "", conda_base)
            conda_exe <- file.path(conda_base, "bin", "conda")
            if (file.exists(conda_exe)) {
              message("Found conda via shell config: ", conda_base)
              return(list(
                conda_base = conda_base,
                conda_exe = conda_exe,
                mamba_exe = file.path(conda_base, "bin", "mamba")
              ))
            }
          }
        }
      }
    }
  }

  stop("Could not find conda installation. Please set CONDA_HOME environment variable or ensure conda is properly installed.")
}

# ---------- Package installation ----------
install_conda_packages <- function(env_name) {
  sys_name <- Sys.info()[["sysname"]]
  is_windows <- sys_name == "Windows"
  is_macos   <- sys_name == "Darwin"

  conda_info <- find_conda_path()

  has_nvidia <- tryCatch({
    system("nvidia-smi", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
  }, error = function(e) FALSE)

  # PyTorch + scvi-tools logic
  if (is_macos) {
    torch_cmd <- "pip3 install torch torchvision"
    scvi_cmd  <- "pip3 install -U \"scvi-tools[metal,scanpy,autotune,interpretability,regseq]\""
  } else if (is_windows) {
    if (has_nvidia) {
      torch_cmd <- "pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu128"
    } else {
      torch_cmd <- "pip3 install torch torchvision"
    }
    scvi_cmd <- if (has_nvidia) {
      "pip3 install -U \"scvi-tools[cuda,scanpy,autotune,interpretability,regseq]\""
    } else {
      "pip3 install -U \"scvi-tools[scanpy,autotune,interpretability,regseq]\""
    }
  } else { # Linux
    if (has_nvidia) {
      torch_cmd <- "pip3 install torch torchvision"
      scvi_cmd  <- "pip3 install -U \"scvi-tools[cuda,scanpy,autotune,interpretability,regseq]\""
    } else {
      torch_cmd <- "pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cpu"
      scvi_cmd  <- "pip3 install -U \"scvi-tools[scanpy,autotune,interpretability,regseq]\""
    }
  }

  # Determine mamba executable path
  if (is_windows) {
    if (is.null(conda_info$mamba_exe) || !file.exists(conda_info$mamba_exe)) {
      warning("Mamba not found, falling back to conda for package management")
      mamba_exe <- "conda"
    } else {
      mamba_exe <- conda_info$mamba_exe
    }
  } else {
    mamba_exe <- conda_info$mamba_exe
    if (!file.exists(mamba_exe)) {
      warning("Mamba not found, falling back to conda for package management")
      mamba_exe <- conda_info$conda_exe
    }
  }

  # Build command - simplified and consistent across platforms
  bash_cmd <- sprintf(
    "bash -c 'source %s/etc/profile.d/conda.sh && \
    conda activate %s && \
    %s install numpy==1.26.4 -y -v && \
    %s && \
    %s && \
    pip3 install jax==0.6.2 cellbender==0.3.0 scrublet==0.2.3 umap-learn==0.5.8 && \
    pip3 install --no-cache-dir -U git+https://github.com/broadinstitute/CellBender.git@4334e8966217c3591bf7c545f31ab979cdc6590d'",
    conda_info$conda_base,
    env_name,
    mamba_exe,
    torch_cmd,
    scvi_cmd
  )

  result <- system(bash_cmd)
  if (result != 0) {
    warning("Package installation may have failed. Exit code: ", result)
  }
  return(result)
}

#' Setup Conda Environment for Single-Cell Analysis
#'
#' @description
#' Internal helper function that creates a new conda environment with specified
#' Python version and installs all required packages for single-cell RNA-seq
#' analysis. Called automatically by \code{\link{prepare.python}} when the
#' environment doesn't exist.
#'
#' @param env_name Character string. Name for the new conda environment
#' @param python_version Character string. Python version to install. Default: "3.12"
#'
#' @return Returns the exit code from the installation command (0 = success)
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Locates conda/mamba installation
#'   \item Creates a new conda environment with specified Python version
#'   \item Calls \code{install_conda_packages()} to install all required packages
#' }
#'
#' @keywords internal
#' @export
setup_conda_env <- function(env_name = "scrna_analysis", python_version = "3.12") {
  sys_name <- Sys.info()[["sysname"]]
  is_windows <- sys_name == "Windows"

  conda_info <- find_conda_path()

  # Determine mamba executable path
  if (is_windows) {
    if (is.null(conda_info$mamba_exe) || !file.exists(conda_info$mamba_exe)) {
      warning("Mamba not found, falling back to conda for package management")
      mamba_exe <- "conda"
    } else {
      mamba_exe <- conda_info$mamba_exe
    }
  } else {
    mamba_exe <- conda_info$mamba_exe
    if (!file.exists(mamba_exe)) {
      warning("Mamba not found, falling back to conda for package management")
      mamba_exe <- conda_info$conda_exe
    }
  }

  # Use -n (name) for all platforms to ensure environment is discoverable by name
  bash_cmd <- sprintf(
    "bash -c 'source %s/etc/profile.d/conda.sh && \
     %s create -n %s python=%s -y'",
    conda_info$conda_base, mamba_exe, env_name, python_version
  )

  result <- system(bash_cmd)
  if (result != 0) {
    stop("Environment creation failed. Exit code: ", result)
  }

  install_conda_packages(env_name)
}

#' Prepare Python Environment for Single-Cell Analysis
#'
#' @description
#' Checks for and optionally creates a conda environment with Python packages
#' required for single-cell RNA-seq analysis (CellBender, scVI, Scanpy, etc.).
#' Automatically detects the operating system and GPU availability to install
#' appropriate package versions. Validates that all required modules are
#' importable and offers to reinstall if modules are missing.
#'
#' @param conda_env Character string. Name of the conda environment to use or
#'   create. Default: "scrna_analysis"
#' @param python_version Character string. Python version to use when creating
#'   a new environment. Default: "3.12"
#'
#' @return Invisibly returns TRUE. Called primarily for its side effects of
#'   environment setup and validation.
#'
#' @details
#' **Environment Detection:**
#'
#' The function automatically detects:
#' \itemize{
#'   \item \strong{Operating system}: Windows, macOS, or Linux
#'   \item \strong{GPU availability}: Checks for NVIDIA GPU using nvidia-smi
#'   \item \strong{Conda installation}: Searches common installation paths and
#'         shell configuration files
#'   \item \strong{Mamba availability}: Prefers mamba for faster package resolution,
#'         falls back to conda if unavailable
#' }
#'
#' **Installed Packages:**
#'
#' The environment includes:
#' \itemize{
#'   \item \strong{Core packages}:
#'     \itemize{
#'       \item NumPy 1.26.4 (pinned for compatibility)
#'       \item PyTorch and torchvision (CPU or CUDA version based on GPU availability)
#'     }
#'   \item \strong{Single-cell analysis}:
#'     \itemize{
#'       \item scvi-tools (with appropriate extras: metal for macOS, cuda for
#'             Linux/Windows with GPU)
#'       \item Scanpy (single-cell analysis)
#'       \item CellBender 0.3.0 + latest development version from GitHub
#'       \item Scrublet 0.2.3 (doublet detection)
#'     }
#'   \item \strong{Supporting packages}:
#'     \itemize{
#'       \item JAX 0.6.2
#'       \item UMAP-learn 0.5.8
#'     }
#' }
#'
#' **Platform-Specific Installation:**
#'
#' \strong{macOS}:
#' \itemize{
#'   \item PyTorch: CPU version (no CUDA support)
#'   \item scvi-tools: metal extras for Apple Silicon GPU acceleration
#' }
#'
#' \strong{Windows}:
#' \itemize{
#'   \item PyTorch: CUDA 12.8 if GPU detected, CPU otherwise
#'   \item scvi-tools: cuda extras if GPU detected
#'   \item Requires CONDA_HOME environment variable set
#' }
#'
#' \strong{Linux}:
#' \itemize{
#'   \item PyTorch: CUDA if GPU detected, CPU otherwise
#'   \item scvi-tools: cuda extras if GPU detected
#' }
#'
#' **Conda Discovery:**
#'
#' The function searches for conda in the following order:
#' \enumerate{
#'   \item \strong{Windows}: CONDA_HOME environment variable (required)
#'   \item \strong{Unix}: Common installation paths:
#'     \itemize{
#'       \item ~/miniconda3, ~/anaconda3
#'       \item /opt/conda, /opt/miniconda3, /opt/anaconda3
#'       \item /usr/local/miniconda3, /usr/local/anaconda3
#'     }
#'   \item \strong{Unix}: Shell configuration files (.bashrc, .bash_profile,
#'         .zshrc, .profile) for conda.sh sourcing
#' }
#'
#' **Interactive Features:**
#'
#' If required modules are missing from an existing environment, the function:
#' \enumerate{
#'   \item Lists the missing modules
#'   \item Prompts the user to reinstall
#'   \item Reinstalls all packages if user confirms
#' }
#'
#' **Helper Functions:**
#'
#' The main function uses internal helpers (not exported):
#' \itemize{
#'   \item \code{find_conda_path()}: Locates conda installation
#'   \item \code{setup_conda_env()}: Creates new conda environment
#'   \item \code{install_conda_packages()}: Installs all required packages
#' }
#'
#' @section First-Time Setup:
#' On first use, the function will:
#' \enumerate{
#'   \item Detect conda installation
#'   \item Create the specified environment with Python 3.12
#'   \item Install all required packages (may take 10-30 minutes)
#'   \item Configure reticulate to use the new environment
#' }
#'
#' @section Subsequent Uses:
#' On subsequent calls, the function:
#' \enumerate{
#'   \item Activates the existing environment
#'   \item Validates all required modules are importable
#'   \item Offers to reinstall if any modules are missing
#' }
#'
#' @section Troubleshooting:
#' \strong{Windows}:
#' \itemize{
#'   \item Set CONDA_HOME: \code{Sys.setenv(CONDA_HOME = "C:/Users/YourName/miniconda3")}
#'   \item Ensure conda/mamba are in System PATH
#' }
#'
#' \strong{macOS/Linux}:
#' \itemize{
#'   \item Ensure conda is initialized in your shell: \code{conda init bash}
#'   \item Check conda.sh is sourced in ~/.bashrc or ~/.zshrc
#' }
#'
#' \strong{All platforms}:
#' \itemize{
#'   \item Install mamba for faster package resolution: \code{conda install mamba -n base -c conda-forge}
#'   \item Check environment exists: \code{conda env list}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage - creates or checks default environment
#' prepare.python()
#'
#' # Use custom environment name
#' prepare.python(conda_env = "my_scrna_env")
#'
#' # Specify different Python version
#' prepare.python(python_version = "3.11")
#'
#' # For Windows, set CONDA_HOME first
#' Sys.setenv(CONDA_HOME = "C:/Users/YourName/miniconda3")
#' prepare.python()
#'
#' # After setup, Python packages are available via reticulate
#' library(reticulate)
#' sc <- import("scanpy")
#' scvi <- import("scvi")
#' cb <- import("cellbender")
#'
#' # Verify GPU availability for PyTorch (if applicable)
#' torch <- import("torch")
#' torch$cuda$is_available()  # Should return TRUE if CUDA GPU detected
#' }
#'
#' @section System Requirements:
#' \itemize{
#'   \item \strong{Required}: Conda or Miniconda installation
#'   \item \strong{Recommended}: Mamba for faster package installation
#'   \item \strong{Optional}: NVIDIA GPU with CUDA support for GPU acceleration
#'   \item \strong{Disk space}: ~5-10 GB for full environment
#'   \item \strong{Internet}: Required for package downloads
#' }
#'
#' @note
#' \itemize{
#'   \item This function is designed for use with \code{\link{run_scvi_integration}}
#'         and other Python-based single-cell tools
#'   \item Package versions are pinned for stability and compatibility
#'   \item The function uses bash commands and may not work in some restricted
#'         environments
#'   \item Environment creation can take 10-30 minutes depending on internet
#'         speed and system resources
#' }
#'
#' @references
#' CellBender: https://cellbender.readthedocs.io/
#'
#' scvi-tools: https://docs.scvi-tools.org/
#'
#' Scanpy: https://scanpy.readthedocs.io/
#'
#' @importFrom reticulate use_condaenv import
#'
#' @seealso
#' \code{\link{run_scvi_integration}} for scVI-based batch correction
#' \code{\link[reticulate]{use_condaenv}} for conda environment management
#' \code{\link[reticulate]{import}} for importing Python modules
#'
#' @export
prepare.python <- function(conda_env = "scrna_analysis", python_version = "3.12") {

  message("🔍 Checking environment '", conda_env, "'...")
  require(reticulate)

  # Configure reticulate to find conda properly
  conda_info <- find_conda_path()
  Sys.setenv(RETICULATE_CONDA = conda_info$conda_exe)

  env_exists <- FALSE

  # Try to use the environment by name (consistent across platforms)
  try({
    use_condaenv(conda_env, required = TRUE)
    env_exists <- TRUE
  }, silent = TRUE)

  if (!env_exists) {
    message("❌ Environment not found. Creating...")
    setup_conda_env(conda_env, python_version)

    # Try to use the environment again after creation
    try({
      use_condaenv(conda_env, required = TRUE)
    }, silent = TRUE)

    return(invisible(TRUE))
  }

  message("✅ Environment found. Testing modules...")
  modules <- c("cellbender", "scanpy", "scrublet", "scvi", "umap")
  missing <- c()

  for (mod in modules) {
    ok <- tryCatch({
      import(mod)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) missing <- c(missing, mod)
  }

  if (length(missing) > 0) {
    message("⚠ Missing modules detected: ", paste(missing, collapse = ", "))
    ans <- readline("Do you want to reinstall missing modules? (y/n): ")
    if (tolower(ans) == "y") {
      message("🔄 Reinstalling missing modules...")
      install_conda_packages(conda_env)
    } else {
      message("ℹ Skipped module reinstallation.")
    }
  } else {
    message("✅ All required modules are present.")
  }

  return(invisible(TRUE))
}
