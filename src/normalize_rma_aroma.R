#' Normalizes big affymetrix experiment using aroma.affymetrix
#' 
#' @details 
#' The `aroma` approach is quite different from `affy`, thus this script is left
#' for manual execution
#' 
#' @references 
#'   http://www.aroma-project.org/howtos/ImproveProcessingTime/
#'   http://www.aroma-project.org/archive/memoization/
#' 
normalize_rma_aroma <- function(raw_path, exp_path, exp_name, chipdb){
  
  # Try to speed up (may want to put into ./.Rprofile file)
  future::plan("sequential") # allow parallel
  setOption(aromaSettings, "memory/ram", 200) # 1000-fold more mem than default
  
  # Clean possibly misleading cache from previous runs
  # Note that this may slow down processing, so could be removed once everything works smoothly
  R.cache::clearCache(recursive = T, prompt = F)
  
  # Prepare directory structures, copy and move files, pull out binary CDF
  chipdb <- .prepare_aroma_structure(exp_path, raw_path, exp_name, chipdb)

  # Prepare aroma objects
  annodb <- aroma.affymetrix::AffymetrixCdfFile$byChipType(chipdb)
  affy_raw <- aroma.affymetrix::AffymetrixCelSet$byName(exp_name, chipType = chipdb)
  
  # do rma and extract expression set
  affy_norm_temp <- aroma.affymetrix::doRMA(affy_raw, flavor = "oligo", uniquePlm = FALSE,  verbose = TRUE)
  affy_norm <- aroma.affymetrix::extractExpressionSet(affy_norm_temp, annotationPkg = NULL)
  
  if (length(affy_norm@annotation) == 0){
    affy_norm@annotation <- chipdb
  }

  # str(affy_norm)
  # head((affy_norm@assayData$exprs))
  # print(paste0("Nubmer of NA in exprs: ", sum(is.na(affy_norm@assayData$exprs))))
  # save("affy_norm", file = "affy_norm.RData") # may want to compress
  
  return(affy_norm)
}

.prepare_aroma_structure <- function(exp_path, raw_path, exp_name, annodb){
  
  chiptype <- tolower(gsub("_|-|cdf$", "", tools::file_path_sans_ext(basename(annodb))))
  print(paste("Chip type:", chiptype))
  print(paste("ChipDB:", annodb))
  
  # Experiment direcotry--------------------------------------------------------
  if (!dir.exists(file.path(exp_path))) {
    print(paste("Creating dir at:", exp_path, "..."))
    dir.create(file.path(exp_path))
  }
  
  setwd(exp_path)
  
  # rawData directory-----------------------------------------------------------
  rawdata_path <- file.path(exp_path,"rawData", exp_name, chiptype)
  if (!dir.exists(file.path(rawdata_path))) {
    print(paste("Creating dir at:", rawdata_path, "..."))
    dir.create(file.path(rawdata_path), recursive = T)
  }
  
  source_files <- list.files(raw_path, "\\.CEL(\\.gz| ?)$", full.names = T)
  target_files <- list.files(rawdata_path, "\\.CEL(\\.gz| ?)$",  full.names = T)
  to_be_copied <- source_files[!(gsub("\\.gz$", "", basename(source_files)) %in%
                                 gsub("\\.gz$", "", basename(target_files)))]
  
  print("Setting up rawData folder (includes moving and gunzipping lot of files), be patient ...")
  if (!(path.expand(raw_path) == path.expand(rawdata_path))){
    # system(paste0("cp -n ", raw_path, "/* ", rawdata_path))
    if (length(to_be_copied) > 0){
      is_ok <- file.copy(to_be_copied, rawdata_path, overwrite = F)
      failed <- to_be_copied[!is_ok]
      if (length(failed) > 0){
        print(paste0("Could not copy ", length(failed), " files (",
                     paste0(failed, collapse = ","), ")."))
      }
    }
  }
  system(paste0("gunzip ", rawdata_path, "/*.gz"))
  
  #annotationData directory-----------------------------------------------------
  print("Setting up annotationData folder ...")
  annodata_path <- file.path(exp_path,"annotationData", "chipTypes", chiptype)
  if (!dir.exists(file.path(annodata_path))) {
    print(paste("Creating dir at:", annodata_path, "..."))
    dir.create(file.path(annodata_path), recursive = T)
  }
  ## creation of annotation library---------------------------------------------
  annodb <- .get_cdf(annodb, annodata_path, rawdata_path, chiptype)
  
  return(annodb)
  
}

#' Get a valid CDF annotation file
#' 
#' @details 
#' Returns name of a valid annotation file and puts it into location expected by aroma
.get_cdf <- function(annodb, annodata_path, rawdata_path, chiptype){
  if (file.exists(annodb)){
    annodb <- .make_binary_cdf_from_file(annodb, annodata_path, chiptype)
  } else {
    annodb <- .get_cdf_from_package(annodb, annodata_path, rawdata_path)
  }
  return(annodb)
}

#' Get a valid CDF annotation file from existing R package
#' 
#' @details 
#' package must be installed (e.g. from bioconductor)
.get_cdf_from_package <- function(annodb, annodata_path, rawdata_path){
  suppressPackageStartupMessages(library(annodb, character.only = T))
  fnames <- list.files(rawdata_path, full.names = T)
  targetf <- file.path(annodata_path, paste0(gsub("cdf$", "",annodb), ".cdf"))
  
  if (!(file.exists(targetf))){
    print("Pulling out CDF from package....")
    pathname <- aroma.affymetrix::env2Cdf(annodb, fnames[[1]])
    system(paste0("mv -v ", pathname, " ", annodata_path))
  } else{
    print(paste0("Using existing cdf at: ",targetf," ...."))
    pathname <- targetf
  }
  
  return(tools::file_path_sans_ext(basename(pathname)))
}

#' Get a valid CDF annotation file from ASCII cdf file 
#' 
#' @details 
#' File can be obtain e.g. from BrainArray website
.make_binary_cdf_from_file <- function(annodb, annodata_path, chiptype){
  
  assertthat::assert_that(file.exists(annodb))
  
  out_path <- file.path(annodata_path, paste0(chiptype, ".cdf"))
  if (file.exists(out_path)){
    print(paste("Using file at:", out_path))
  } else {
    print("Converting ASCII CDF to binary....")
    affxparser::convertCdf(annodb, out_path, force = TRUE)
  }
  
  return(tools::file_path_sans_ext(basename(out_path)))
}

#' Remove CEL files from experiment location
#' 
#' @details 
#' Does not remove if the raw files are present only once
#' 
.cleanup_rawdata <- function(rawdata_path, raw_path){
  # cleans up a rawData folder in experiment location
  celfiles <- affy::list.celfiles(rawdata_path, full.names = T)
  if (!(path.expand(raw_path) == path.expand(rawdata_path))){
    print(paste0("Removing CEL files from: ", rawdata_path))
    file.remove(celfiles)
  }
  return(invisible(TRUE))
}