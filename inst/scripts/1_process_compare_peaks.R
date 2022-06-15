start <- Sys.time()

library(package = baseline, quietly = TRUE)
library(package = data.table, quietly = TRUE)
library(package = docopt, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = future, quietly = TRUE)
library(package = future.apply, quietly = TRUE)
library(package = MSnbase, quietly = TRUE)
library(package = nucleR, quietly = TRUE)
library(package = parallel, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = pracma, quietly = TRUE)
library(package = progressr, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = xcms, quietly = TRUE)

source(file = "R/change_intensity_name.R")
source(file = "R/check_export_dir.R")
source(file = "R/compare_peaks.R")
source(file = "R/extract_ms.R")
source(file = "R/extract_ms_peak.R")
source(file = "R/filter_ms.R")
source(file = "R/get_gnps.R")
source(file = "R/get_params.R")
source(file = "R/improve_signal.R")
source(file = "R/improve_signals_progress.R")
source(file = "R/join_peaks.R")
source(file = "R/log_debug.R")
source(file = "R/make_confident.R")
source(file = "R/normalize_chromato.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/peaks_progress.R")
source(file = "R/plot_chromatogram.R")
source(file = "R/prepare_features.R")
source(file = "R/prepare_mz.R")
source(file = "R/prepare_peaks.R")
source(file = "R/prepare_rt.R")
source(file = "R/preprocess_chromatograms.R")
source(file = "R/preprocess_peaks.R")
source(file = "R/process_peaks.R")
source(file = "R/transform_baseline.R")
source(file = "R/transform_ms.R")
source(file = "R/y_as_na.R")

progressr::handlers(global = TRUE)
progressr::handlers("progress")

step <- "processing"
paths <- parse_yaml_paths()
params <- ""
params <- get_params(step = step)

log_debug(
  "This program performs",
  "Quantitative and Qualitative Contextualization",
  "of in depth annotated extracts"
)
log_debug("Authors: \n", "AR")
log_debug("Contributors: \n", "...")

#' Dirty generic paths and parameters
source(file = "R/dirty_paths_params.R")
#' Specific paths

log_debug(x = "listing files")
files <- list.files(
  path = TOYSET,
  pattern = ".mzML",
  full.names = TRUE,
  recursive = TRUE
)

# files <- files[grepl(pattern = "M_17|M_28|M_36|M_40|M_47|M_57|M_67", x = files)]

names <- list.files(
  path = TOYSET,
  pattern = ".mzML",
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML",
    replacement = "",
    fixed = TRUE
  )

# names <- names[grepl(pattern = "M_17|M_28|M_36|M_40|M_47|M_57|M_67", x = names)]

log_debug(x = "loading feature table")
# feature_table <- read_features(id = GNPS_JOB)
feature_table <- readr::read_delim(file = FEATURES)

log_debug(x = "loading raw files (can take long if loading multiple files)")
dda_data <- MSnbase::readMSData(
  files = files,
  mode = "onDisk",
  msLevel. = 1
)

log_debug(x = "opening raw files objects and extracting chromatograms")
chromatograms_all <- lapply(files, mzR::openMSfile) |>
  lapply(mzR::chromatograms) |>
  purrr::flatten()

log_debug(x = "preparing feature list ...")
df_features <- feature_table |>
  prepare_features()

if (params$signal$detector$bpi == TRUE) {
  detector <- "bpi"
}
if (params$signal$detector$bpi == TRUE) {
  chromatograms_list_bpi <- preprocess_chromatograms(detector = "bpi")
  peaks_prelist_bpi <- preprocess_peaks(detector = "bpi")
  peaks_list_pda <- process_peaks(detector = "bpi")
}
if (params$signal$detector$cad == TRUE) {
  detector <- "cad"
}
if (params$signal$detector$cad == TRUE) {
  chromatograms_list_cad <- preprocess_chromatograms()
  peaks_prelist_cad <- preprocess_peaks()
  peaks_list_cad <- process_peaks()
}
if (params$signal$detector$pda == TRUE) {
  detector <- "pda"
}
if (params$signal$detector$pda == TRUE) {
  #' TODO check if change for PDA doing baselining on normal and not improved
  chromatograms_list_pda <-
    preprocess_chromatograms(detector = "pda")
  peaks_prelist_pda <- preprocess_peaks(detector = "pda")
  peaks_list_pda <- process_peaks(detector = "pda")
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
