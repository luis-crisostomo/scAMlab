# Function 1: Prepare data for CellBender (run locally) - Modified to output .h5
#' @export
PrepareCellBender_remote <- function(folder.path, expected.cells = 5000, total.droplets = 25000,
                                     learning.rate = 0.0001, fpr = 0.01, epochs = 150, 
                                     remote.path = NULL, email = NULL) {
  
  message("Preparing data for CellBender processing...")
  
  # Validate input path
  if (!dir.exists(folder.path)) {
    stop(paste("Input directory does not exist:", folder.path))
  }
  
  original_path <- normalizePath(folder.path, mustWork = TRUE)
  sample_name <- basename(dirname(original_path))
  
  message(paste("Processing sample:", sample_name))
  message(paste("Input directory:", original_path))
  
  # Create output directory for H5 file
  parent_dir <- dirname(original_path)
  output_h5_dir <- file.path(parent_dir, "for_cellbender_h5")
  dir.create(output_h5_dir, showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    
    # STEP 1: Find Parse files
    message("Step 1: Finding Parse Biosciences files...")
    all_files <- list.files(original_path, full.names = TRUE, recursive = TRUE)
    
    # Find the three required files
    cell_metadata_file <- all_files[grepl("cell_metadata|barcodes", basename(all_files), ignore.case = TRUE)][1]
    gene_metadata_file <- all_files[grepl("all_genes|features|genes", basename(all_files), ignore.case = TRUE)][1]
    count_matrix_file <- all_files[grepl("count_matrix|matrix", basename(all_files), ignore.case = TRUE)][1]
    
    if (is.na(cell_metadata_file) || is.na(gene_metadata_file) || is.na(count_matrix_file)) {
      stop("Could not find required Parse files (cell_metadata, all_genes, count_matrix)")
    }
    
    message(paste("Found files:"))
    message(paste("  Cell metadata:", basename(cell_metadata_file)))
    message(paste("  Gene metadata:", basename(gene_metadata_file)))
    message(paste("  Count matrix:", basename(count_matrix_file)))
    
    # STEP 2: Handle gzipped files
    message("Step 2: Checking for gzipped files...")
    
    temp_base <- tempdir()
    temp_parse_dir <- file.path(temp_base, paste0(sample_name, "_parse_temp"))
    dir.create(temp_parse_dir, showWarnings = FALSE, recursive = TRUE)
    
    files_to_use <- list(
      cell_metadata = cell_metadata_file,
      gene_metadata = gene_metadata_file,
      count_matrix = count_matrix_file
    )
    
    target_names <- c(
      cell_metadata = "cell_metadata.csv",
      gene_metadata = "all_genes.csv",
      count_matrix = "count_matrix.mtx"
    )
    
    use_temp_dir <- FALSE
    
    for (file_type in names(files_to_use)) {
      source_file <- files_to_use[[file_type]]
      
      if (grepl("\\.gz$", source_file)) {
        message(paste("Extracting gzipped file:", basename(source_file)))
        target_file <- file.path(temp_parse_dir, target_names[[file_type]])
        R.utils::gunzip(source_file, destname = target_file, remove = FALSE, overwrite = TRUE)
        files_to_use[[file_type]] <- target_file
        use_temp_dir <- TRUE
      }
    }
    
    # Determine which directory to use for ReadParseBio
    parse_dir <- if (use_temp_dir) {
      # Copy non-gzipped files to temp directory if some files were extracted
      for (file_type in names(files_to_use)) {
        source_file <- files_to_use[[file_type]]
        if (!grepl(temp_parse_dir, source_file)) {
          target_file <- file.path(temp_parse_dir, target_names[[file_type]])
          file.copy(source_file, target_file, overwrite = TRUE)
        }
      }
      temp_parse_dir
    } else {
      original_path
    }
    
    # STEP 3: Assemble Seurat object with ReadParseBio
    message("Step 3: Reading Parse data with ReadParseBio...")
    counts_data <- Seurat::ReadParseBio(data.dir = parse_dir)
    message(paste("Loaded matrix:", nrow(counts_data), "genes x", ncol(counts_data), "cells"))
    
    # STEP 4: Create H5 file for CellBender
    message("Step 4: Creating H5 file for CellBender...")
    
    # Define output H5 file path
    h5_output_file <- file.path(output_h5_dir, paste0(sample_name, "_raw_feature_bc_matrix.h5"))
    
    # Write to H5 format using DropletUtils
    DropletUtils::write10xCounts(h5_output_file, counts_data,
                                 barcodes = colnames(counts_data),
                                 gene.id = rownames(counts_data),
                                 gene.symbol = rownames(counts_data),
                                 type = "HDF5",  # H5 format
                                 overwrite = TRUE)
    
    # Verify H5 file was created
    if (file.exists(h5_output_file)) {
      message(paste("Successfully created H5 file:", basename(h5_output_file)))
      message(paste("File size:", round(file.size(h5_output_file) / 1024^2, 2), "MB"))
    } else {
      stop("Failed to create H5 file")
    }
    
    # Test reading the H5 file to ensure it's valid
    message("Step 5: Validating H5 file...")
    tryCatch({
      test_read <- Seurat::Read10X_h5(h5_output_file)
      message(paste("H5 file validation successful:", nrow(test_read), "genes x", ncol(test_read), "cells"))
    }, error = function(e) {
      stop(paste("H5 file validation failed:", e$message))
    })
    
    # Create metadata file for remote processing
    metadata <- list(
      sample_name = sample_name,
      original_cells = ncol(counts_data),
      original_genes = nrow(counts_data),
      expected_cells = expected.cells,
      total_droplets = total.droplets,
      learning_rate = learning.rate,
      fpr = fpr,
      epochs = epochs,
      input_path = original_path,
      output_h5_dir = output_h5_dir,
      remote_path = if(is.null(remote.path)) "" else remote.path,
      h5_file = h5_output_file,
      processing_date = Sys.time()
    )
    
    saveRDS(metadata, file.path(output_h5_dir, "metadata.rds"))
    
    # Create a comprehensive shell script for CellBender execution
    cellbender_script <- file.path(output_h5_dir, paste0("run_cb_", sample_name, ".sh"))
    
    # Build email notification lines
    email_lines <- if (!is.null(email)) {
      paste0("#SBATCH --mail-type=END,FAIL\n",
             "#SBATCH --mail-user=", email, "\n")
    } else {
      ""
    }
    
    # Set paths - UPDATED: changed "cellbender_filtered" to "postCB"
    remote_base <- if(is.null(remote.path) || remote.path == "") "/scratch/project_2007503" else remote.path
    input_path <- paste0(remote_base, "/", sample_name, "/", basename(h5_output_file))
    output_dir <- paste0(remote_base, "/", sample_name, "_cbOutput")
    output_file <- paste0(output_dir, "/", sample_name, "_postCB.h5")
    
    script_content <- paste0(
      "#!/bin/bash\n",
      "#SBATCH --account=project_2007503\n",
      "#SBATCH --partition=gpu\n",
      "#SBATCH --time=02:00:00\n",
      "#SBATCH --mem=32G\n",
      "#SBATCH --gres=gpu:v100:1\n",
      "#SBATCH --cpus-per-task=4\n",
      email_lines,
      "\n",
      "# Print job info\n",
      "echo \"Job started at: $(date)\"\n",
      "echo \"Job ID: $SLURM_JOB_ID\"\n",
      "echo \"Node: $SLURM_NODELIST\"\n",
      "echo \"Working directory: $(pwd)\"\n",
      "echo \"Sample: ", sample_name, "\"\n\n",
      
      "# Set working directory in the same partition as the output\n",
      "cd \"/scratch/project_2007503\"\n\n",
      
      "# Set temporary directory with unique job ID\n",
      "export TMPDIR=/scratch/project_2007503/cellbender_tmp_${SLURM_JOB_ID}\n",
      "mkdir -p $TMPDIR\n",
      "echo \"Temporary directory created: $TMPDIR\"\n\n",
      
      "# Create and set up output directory\n",
      "OUTPUT_DIR=\"", output_dir, "\"\n",
      "mkdir -p $OUTPUT_DIR\n",
      "chmod 755 $OUTPUT_DIR\n",
      "echo \"Output directory created: $OUTPUT_DIR\"\n\n",
      
      "# Load modules and activate environment\n",
      "echo \"Loading modules and activating environment...\"\n",
      "module load bioconda\n",
      "eval \"$(conda shell.bash hook)\"\n",
      "conda activate cellbender\n",
      "module load cuda/11.7.0\n\n",
      #"module load pytables\n\n",  # Load pytables module for ptrepack
      
      "# Verify cellbender is available\n",
      "echo \"Checking CellBender installation...\"\n",
      "which cellbender\n",
      "cellbender --help > /dev/null 2>&1\n",
      "if [ $? -ne 0 ]; then\n",
      "    echo \"ERROR: CellBender is not properly installed or accessible\"\n",
      "    exit 1\n",
      "fi\n\n",
      
      "# Verify Python modules are available\n",
      "python --version\n",
      "python -c \"import tables; print('PyTables:', tables.__version__)\"\n",
      "python -c \"import torch; print('CUDA available:', torch.cuda.is_available())\"\n\n",
      
      "# Verify input file exists\n",
      "INPUT_FILE=\"", input_path, "\"\n",
      "if [ ! -f \"$INPUT_FILE\" ]; then\n",
      "    echo \"ERROR: Input file $INPUT_FILE not found\"\n",
      "    echo \"Current directory contents:\"\n",
      "    ls -la\n",
      "    exit 1\n",
      "fi\n",
      "echo \"Input file found: $INPUT_FILE\"\n\n",
      
      "echo \"Starting CellBender at: $(date)\"\n",
      "echo \"Parameters:\"\n",
      "echo \"  Expected cells: ", expected.cells, "\"\n",
      "echo \"  Total droplets: ", total.droplets, "\"\n",
      "echo \"  Learning rate: ", learning.rate, "\"\n",
      "echo \"  FPR: ", fpr, "\"\n",
      "echo \"  Epochs: ", epochs, "\"\n\n",
      
      "# UPDATED: Set checkpoint directory to temporary folder\n",
      "export CELLBENDER_CHECKPOINT_DIR=$TMPDIR\n",
      "echo \"CellBender checkpoint directory set to: $CELLBENDER_CHECKPOINT_DIR\"\n\n",
      
      "# Run CellBender\n",
      "cellbender remove-background \\\n",
      "    --cuda \\\n",
      "    --input \"$INPUT_FILE\" \\\n",
      "    --output \"", output_file, "\" \\\n",
      "    --expected-cells ", expected.cells, " \\\n",
      "    --total-droplets-included ", total.droplets, " \\\n",
      "    --learning-rate ", learning.rate, " \\\n",
      "    --fpr ", fpr, " \\\n",
      "    --epochs ", epochs, "\n\n",
      
      "CELLBENDER_EXIT_CODE=$?\n",
      "echo \"CellBender finished at: $(date)\"\n",
      "echo \"CellBender exit code: $CELLBENDER_EXIT_CODE\"\n\n",
      
      "# Check if output was created\n",
      "FILTERED_FILE=\"", paste0(output_dir, "/", sample_name, "_postCB_filtered.h5"), "\"\n",
      "if [ -f \"$FILTERED_FILE\" ]; then\n",
      "    echo \"SUCCESS: CellBender output file created successfully\"\n",
      "    ls -lh \"$OUTPUT_DIR\"/\n",
      "else\n",
      "    echo \"ERROR: CellBender output file was not created\"\n",
      "    exit 1\n",
      "fi\n\n",
      
      "# UPDATED: Convert H5 file to Seurat-compatible format using ptrepack\n",
      "echo \"Converting H5 file to Seurat-compatible format...\"\n",
      "SEURAT_FILE=\"", paste0(output_dir, "/", sample_name, "_postCB_seurat.h5"), "\"\n",
      "ptrepack --complevel 5 \"$FILTERED_FILE:/matrix\" \"$SEURAT_FILE:/matrix\" --overwrite-nodes\n\n",
      
      "PTREPACK_EXIT_CODE=$?\n",
      "echo \"ptrepack exit code: $PTREPACK_EXIT_CODE\"\n\n",
      
      "if [ $PTREPACK_EXIT_CODE -eq 0 ] && [ -f \"$SEURAT_FILE\" ]; then\n",
      "    echo \"SUCCESS: Seurat-compatible H5 file created successfully\"\n",
      "    ls -lh \"$SEURAT_FILE\"\n",
      "else\n",
      "    echo \"ERROR: Failed to create Seurat-compatible H5 file\"\n",
      "    exit 1\n",
      "fi\n\n",
      
      "# Final output summary\n",
      "echo \"Final output files:\"\n",
      "echo \"  Original CellBender output: $FILTERED_FILE\"\n",
      "echo \"  Seurat-compatible file: $SEURAT_FILE\"\n",
      "ls -lh \"$OUTPUT_DIR\"/\n\n",
      
      "# Clean up temporary files\n",
      "echo \"Cleaning up temporary directory: $TMPDIR\"\n",
      "rm -rf $TMPDIR\n\n",
      
      "echo \"Job completed successfully at: $(date)\"\n"
    )
    
    writeLines(script_content, cellbender_script)
    Sys.chmod(cellbender_script, mode = "0755")  # Make executable
    
    message(paste("Success! H5 file prepared for CellBender in:", output_h5_dir))
    message(paste("CellBender script created:", basename(cellbender_script)))
    message("Key improvements:")
    message("  - Output files will be named with 'postCB' instead of 'cellbender_filtered'")
    message("  - Checkpoint files will be stored in temporary directory")
    message("  - Seurat-compatible H5 file will be created automatically")
    if (!is.null(email)) {
      message(paste("Email notifications will be sent to:", email))
    }
    message("Upload this directory to Puhti for CellBender processing.")
    
    return(list(
      output_dir = output_h5_dir,
      h5_file = h5_output_file,
      metadata = metadata,
      cellbender_script = cellbender_script
    ))
    
  }, finally = {
    # Clean up temporary directories
    if (exists("temp_parse_dir") && dir.exists(temp_parse_dir)) {
      unlink(temp_parse_dir, recursive = TRUE)
    }
  })
}
