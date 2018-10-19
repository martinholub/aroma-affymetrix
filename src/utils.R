.firstup <- function(x) {
  # Make first letter upercase
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# .get_platform <- function(pname = "agilent"){
#   # not used
#   if (startsWith(tolower(pname), "a")){
#     platforms <- list(hs = "hgug4112a.db", mm = "mgug4122a.db", rn = "rgug4131a.db")
#   }else if (startsWith(tolower(pname), "i")){
#     platforms <- list(hs = "lumiHumanIDMapping", mm = "lumiMouseIDMapping",
#                       rn = "lumiRatIDMapping")
#   } else {
#     warning(paste0("get_annotation_lib:wrong-platform: ",
#                    "unknown platform type, selecting agilent"))
#     platforms <- list(hs = "hgug4112a.db", mm = "mgug4122a.db", rn = "rgug4131a.db")
#   }
# }

# .drop_bad_stuf <- function(norm, raw, targets){
#   # This is not used and only example of convenience function
#   bad_samps <- c("GSM738509_RPEChoroid_77_raw.txt")
#   norm <- norm[, colnames(norm) != bad_samps]
#   raw <- raw[colnames(raw) != bad_samps]
#   targets <- targets[rownames(targets) != bad_samps, ]
# }

.pollute_data <- function(data){
  # Polute data with gaussian noise, with exisitng mean and proportional variance
  means <- colMeans(data)
  sds <- sqrt(means) * 0.1
  noise <- matrix(rnorm(prod(dim(data)), mean = rep(0, length(means)), sd = sds), ncol = dim(data)[2])
  data <- data + noise
  return(data)
}

.quick_boxplot <- function(abatch){
  dev.new(noRStudioGD = T)
  # make quick helper plot
  boxplot(log2(data.frame(abatch@assayData$exprs)), outline = F)
  readline(prompt="Press [enter] to continue")
  graphics.off()
}

# Clean Memory in loop
.clean_memory <- function(n=5) { for (i in 1:n) gc(verbose = F) }

#' Fetch config as yaml file 
#' 
#' @description
#' Obtains config either from experiment directory (top-level) or from a default
#' location in `<scripts_dir>/configs/config.yml`. In experiment diretory,
#' both `<exp_name>_config.yml` and `config.yml` are valid names, given that 
#' `<exp_name>` is not NULL.
#' 
#' @param cwd experiment location, also current working directory
#' @param script_dir location of maQCN_pipeline.
#' @param exp_name name of experiment
#' 
#' @return cfg_path path to config
#' @export
#'  
get_config_path <- function(cwd, script_dir, exp_name = NULL){

  if (length(exp_name)>0){
    pattern <- paste0("(",exp_name,"_| ?)config.yml")
  } else {
    pattern <- "config.yml"
  } #   cwd_config <- normalizePath(file.path(cwd, "config.yml"))
  cwd_config <- normalizePath(list.files(file.path(cwd), pattern, full.names = T)[1],
                              mustWork = FALSE)
  
  default_config <- normalizePath(file.path(script_dir, "configs", "config.yml"))
  if (file.exists(cwd_config)) {
    cfg_path <- cwd_config
  } else if (file.exists(default_config)){
    cfg_path <- default_config
  } else {
    stop("No 'config.yml' found in neither experiment nor script location. Aborting.")
  }
  print(paste("maQCNpipe:get_config_path: Using config from", cfg_path, "."))
  
  return(cfg_path)
}

#' Obtain all current warnings and write them to file
#' 
#' @description 
#' Warnings are written to 'WARNINGS.log' file in top-level experiment directory (cwd).
#' If no warnings present, no file is written.
#' 
#' @param location directory where to place the log file, default cwd
#' 
#' @return invisible(TRUE)
#' @export
#' 
make_warnings_file <- function(location = ".", warns = list()){
  if (length(location)==0) {location <- "."}
  warn_file <- normalizePath(file.path(location, "WARNINGS.log"), mustWork = FALSE)
  # warns <- warnings() # assing warnings to variable
  
  if (length(warns)>0){
    print(paste("There were",length(warns),"warnings, see them in WARNINGS.log file."))
    
    con <- file(warn_file, "w")
    sink(con, type = c("output", "message")) # redirect to file
    print(warns)
    sink() # redirect back to stdout,stderr
    close(con) # may not be necessary, but not harmful
  }
  return(invisible(TRUE))
}

#' Strip extensions from character vector of filenames
#' 
#' @details 
#' This slightly duplicates .strip_extensions which is used for single full path
#' TODO: Combine the two once branches merged
#' 
#' @param lin list or vector of characters
#' @return lout character vector
#' @seealso .strip_extensions
.strip_extensions_vec <- function(lin){
  lout <- gsub("(.*)\\.(txt|gpr|cel)( ?|\\.proc)( ?|\\.gz|\\.zip)$", "\\1", lin,
          perl = TRUE, ignore.case = T)
  return(lout)
}