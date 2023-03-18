start <- Sys.time()

library(package = baseline, quietly = TRUE)
library(package = BiocParallel, quietly = TRUE)
library(package = chromatographR, quietly = TRUE)
library(package = data.table, quietly = TRUE)
library(package = docopt, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = future, quietly = TRUE)
library(package = future.apply, quietly = TRUE)
library(package = MSnbase, quietly = TRUE)
library(package = nucleR, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = pracma, quietly = TRUE)
library(package = progressr, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = xcms, quietly = TRUE)

source(file = "R/baseline_chromatogram.R")
source(file = "R/change_intensity_name.R")
source(file = "R/check_export_dir.R")
source(file = "R/compare_peaks.R")
source(file = "R/extract_ms_progress.R")
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
source(file = "R/progressr.R")
source(file = "R/transform_ms.R")
source(file = "R/y_as_na.R")

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
  pattern = paste0(params$filename$mzml, ".mzML"),
  full.names = TRUE,
  recursive = TRUE
)

# files <- files[grepl(pattern = "M_17|M_28|M_36|M_40|M_47|M_57|M_67", x = files)]

names <- list.files(
  path = TOYSET,
  pattern = paste0(params$filename$mzml, ".mzML"),
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
feature_table <- read_delim(file = FEATURES)

log_debug(x = "loading raw files (can take long if loading multiple files)")
dda_data <- readMSData(
  files = files,
  mode = "onDisk",
  msLevel. = 1
)

log_debug(x = "opening raw files objects and extracting chromatograms")
chromatograms_all <- lapply(files, openMSfile) |>
  lapply(chromatograms) |>
  flatten()

log_debug(x = "preparing feature list ...")
df_features <- feature_table |>
  prepare_features()

chromatograms_list_cad <- preprocess_chromatograms()
peaks_prelist_cad <- preprocess_peaks()
detector <- "cad"
peaks_prelist <- switch(detector,
  "bpi" = peaks_prelist_bpi,
  "cad" = peaks_prelist_cad,
  "pda" = peaks_prelist_pda
)
log_debug(x = "processing", detector, "peaks")
log_debug(x = "extracting ms chromatograms (longest step)")
log_debug(x = "count approx 1 minute per 500 features (increasing with features number)")
log_debug(x = "varies a lot depending on features distribution")
list_ms_chromatograms <-
  extract_ms_progress(xs = seq_along(peaks_prelist$list_df_features_with_peaks_long))

log_debug(x = "transforming ms chromatograms")
list_ms_chromatograms_transformed <-
  future_lapply(
    X = list_ms_chromatograms,
    FUN = transform_ms
  )

log_debug(x = "extracting ms peaks")
list_ms_peaks <-
  future_lapply(
    X = list_ms_chromatograms_transformed,
    FUN = extract_ms_peak
  )

log_debug(x = "comparing peaks")
list_comparison_score <-
  future_lapply(
    X = seq_along(list_ms_peaks),
    FUN = compare_peaks
  )

log_debug(x = "selecting features with peaks")
df_features_with_peaks <-
  peaks_prelist$list_df_features_with_peaks_long |>
  bind_rows()

log_debug(x = "There are", nrow(df_features_with_peaks), "features with peaks")

log_debug(x = "summarizing comparison scores")
comparison_scores <- list_comparison_score |>
  flatten()

log_debug(x = "There are", length(comparison_scores), "scores calculated")

log_debug(x = "joining")
df_features_with_peaks$comparison_score <-
  as.numeric(comparison_scores)

log_debug(x = "final aesthetics")
df_features_with_peaks_scored <- df_features_with_peaks |>
  select(
    sample = id,
    peak_id,
    peak_rt_min = rt_min,
    peak_rt_apex = rt_apex,
    peak_rt_max = rt_max,
    peak_area = integral,
    feature_id,
    feature_rt = rt,
    feature_mz = mz,
    feature_area = intensity,
    comparison_score
  ) |>
  distinct()

df_features_without_peaks_scored <-
  peaks_prelist$df_features_without_peaks |>
  mutate(comparison_score = NA) |>
  select(
    sample = id,
    peak_id,
    peak_rt_min = rt_min,
    peak_rt_apex = rt_apex,
    peak_rt_max = rt_max,
    peak_area = integral,
    feature_id,
    feature_rt = rt,
    feature_mz = mz,
    feature_area = intensity,
    comparison_score
  ) |>
  distinct()

log_debug(x = "checking export directory")
check_export_dir(EXPORT_DIR)

log_debug(x = "exporting to ...")
log_debug(x = EXPORT_DIR)
write_tsv(
  x = df_features_with_peaks_scored,
  file = file.path(EXPORT_DIR, switch(detector,
    "bpi" = EXPORT_FILE_BPI,
    "cad" = EXPORT_FILE_CAD,
    "pda" = EXPORT_FILE_PDA
  ))
)
write_tsv(
  x = df_features_without_peaks_scored,
  file = file.path(EXPORT_DIR, switch(detector,
    "bpi" = EXPORT_FILE_BPI_2,
    "cad" = EXPORT_FILE_CAD_2,
    "pda" = EXPORT_FILE_PDA_2
  ))
)

# if (params$signal$detector$pda == TRUE) {}
# TODO check if change for PDA doing baselining on normal and not improved
chromatograms_list_pda <-
  preprocess_chromatograms(
    detector = "pda",
    list = chromatograms_all[c(FALSE, TRUE, FALSE)],
    signal_name = "PDA.1_TotalAbsorbance_0",
    shift = PDA_SHIFT
  )
peaks_prelist_pda <- preprocess_peaks(
  detector = "pda",
  list = chromatograms_list_pda$chromatograms_improved,
  df_long = chromatograms_list_pda$chromatograms_improved_long
)
detector <- "pda"
peaks_list_pda <- process_peaks(detector = "pda")

end <- Sys.time()

log_debug("Script finished in", format(end - start))
