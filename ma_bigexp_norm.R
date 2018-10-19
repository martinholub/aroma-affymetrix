#! /usr/bin/Rscript

#Example call:
# Rscript ma_bigexp_norm.R path/to/config_aroma.yml
# Get Config -------------------------------------------------------------------
cl_args <- commandArgs(trailingOnly = FALSE)
this_dir <- tryCatch(
  {
    this_dir <- dirname(sub("--file=", "", cl_args[grep("--file=", cl_args)]))
  },
    error = function(cond){
      this_dir <- getwd()
  },
  finally = {
      #pass
  })
print(paste0("Script Directory: ", this_dir))
# Retain only arguments passed from command line
cl_args <- cl_args[(grep("--args", cl_args)+1):length(cl_args)]

if (!(length(cl_args) == 0)){
  cfg_path <- cl_args[1]
} else{
  cfg_path <- file.path(this_dir, "configs", "config_aroma.yml")
}
cfg <- config::get(file = cfg_path, use_parent = F)

# Source functions for execution -----------------------------------------------
# TODO: Package this
scripts_dir <- file.path(this_dir, "maQCNpipe", "R")
scripts <- list.files(scripts_dir, pattern = "*.R", full.names = T) 
for (scr in scripts) {source (scr)}
library(aroma.affymetrix)
# Get and Normalize Execution Parameters ---------------------------------------
if ((length(cfg$exp_name) == 0)){
  cfg$exp_name <- basename(cfg$raw_path)
}

if ((length(cfg$chipdb) == 0)){
  # gets cdf package name by looking at first file in the folder
  cfg$chipdb <- .get_annotation_affy(cfg$raw_path)
}

# Execute ----------------------------------------------------------------------
print("Starting aroma.affymetrix big experiment normalization ...")
start_time <- proc.time()
gse_targets <- make_targets(cfg$raw_path)

# Within-experiment normalization (replaces affy::rma)
gse_sum <- normalize_rma_aroma(cfg$raw_path, cfg$exp_path, cfg$exp_name, cfg$chipdb)
rownames(gse_targets) <- colnames(gse_sum)
gse_norm <- build_eset( gse_sum, gse_targets, annodb = gse_sum@annotation,
                        exp_name = cfg$exp_name)

rm(gse_sum)
.clean_memory()
# Compendium normalization, same as before -------------------------------------
norm_files_dir <- file.path(cfg$exp_path, "normData", gse_norm@annotation, "TXT/")
normalize_compendium(gse_norm, norm_files_dir, cfg$rmaTrim, cfg$rmaTv, verbose = TRUE)

# Clean up----------------------------------------------------------------------
# Remove raw data from experiment folder (keeps them in original location)
.cleanup_rawdata(file.path(cfg$exp_path, "rawData", cfg$exp_name, gse_norm@annotation), cfg$raw_path)

## Done routines ---------------------------------------------------------------
make_warnings_file(cfg$exp_path)
end_time <- proc.time() - start_time
print(paste0("Done in: ", end_time[["elapsed"]], "s."))
sessionInfo()