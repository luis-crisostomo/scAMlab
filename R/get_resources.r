#' @export
# Helper function to retrieve system resources on the fly
get_resources <- function() {
  sys <- Sys.info()[["sysname"]]
  
  # --- CPU info ---
  cores <- parallel::detectCores(logical = FALSE)   # physical cores
  logical_cores <- parallel::detectCores(logical = TRUE)  # logical processors
  
  ram <- switch(sys,
                "Linux"   = {
                  out <- system("grep MemTotal /proc/meminfo", intern = TRUE)
                  as.numeric(gsub("[^0-9]", "", out)) / 1024  # in MB
                },
                "Darwin"  = {
                  out <- system("sysctl -n hw.memsize", intern = TRUE)
                  as.numeric(out) / (1024^2)  # in MB
                },
                "Windows" = {
                  out <- system("wmic computersystem get TotalPhysicalMemory", intern = TRUE)
                  bytes <- suppressWarnings(as.numeric(gsub("[^0-9]", "", out)))
                  bytes[!is.na(bytes)] / (1024^2)  # in MB
                },
                NA
  )
  
  cpu_info <- list(
    cpu_cores = cores,
    logical_processors = logical_cores,
    ram_mb = ram
  )
  
  # --- Try NVIDIA GPU ---
  gpu_name <- try(system("nvidia-smi --query-gpu=name --format=csv,noheader", intern = TRUE), silent = TRUE)
  gpu_mem  <- try(system("nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits", intern = TRUE), silent = TRUE)
  
  if (!inherits(gpu_mem, "try-error") && length(gpu_mem) > 0) {
    gpu_info <- list(
      gpu_name = gpu_name,
      gpu_memory_mb = as.numeric(gpu_mem)
    )
    return(c(list(os = sys), cpu_info, gpu_info))
  }
  
  # --- If no NVIDIA GPU ---
  return(c(list(os = sys), cpu_info))
}