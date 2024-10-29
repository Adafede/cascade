start <- Sys.time()

# library(package = chromatographR, quietly = TRUE)
## Getting it this way only to avoid
## "Found more than one class "cluster" in cache issue
source(file = "https://raw.githubusercontent.com/ethanbass/chromatographR/master/R/utils.R")
source(file = "https://raw.githubusercontent.com/ethanbass/chromatographR/master/R/fit_peaks.R")
source(file = "https://raw.githubusercontent.com/ethanbass/chromatographR/master/R/get_peaks.R")
source(file = "https://raw.githubusercontent.com/ethanbass/chromatographR/master/R/get_purity.R")
source(file = "https://raw.githubusercontent.com/ethanbass/chromatographR/master/R/preprocess.R")
source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/log_debug.R")
source(file = "R/check_export_dir.R")
source(file = "R/compare_peaks.R")
source(file = "R/extract_ms_progress.R")
source(file = "R/extract_ms_peak.R")
source(file = "R/prepare_features.R")
source(file = "R/preprocess_chromatograms.R")
source(file = "R/preprocess_peaks.R")
source(file = "R/transform_ms.R")
source(file = "R/zzz.R")

log_debug(
  "This program performs",
  "Quantitative and Qualitative Contextualization",
  "of in depth annotated extracts"
)
log_debug("Authors: \n", "AR")
log_debug("Contributors: \n", "...")

#' Specific paths
AREA_MIN <- 0.005
CAD_SHIFT <- 0.05
INTENSITY_MS_MIN <- 10000
TIME_MIN <- 0.7
TIME_MAX <- 35.2
EXPORT_DIR <- "data/interim/peaks"
EXPORT_FILE_CAD <- "210619_AR_06_V_03_2_01_featuresInformed_cad.tsv.gz"
EXPORT_FILE_CAD_2 <- "210619_AR_06_V_03_2_01_featuresNotInformed_cad.tsv.gz"
FEATURES <- "~/Documents/papers/sapid/sapere_tmp/extract_mzmine/extract.csv"
FILE_POSITIVE <- "data/source/mzml/210619_AR_06_V_03_2_01.mzML"
names <- FILE_POSITIVE |>
  gsub(pattern = ".*/", replacement = "") |>
  gsub(pattern = "[0-9]{8}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML",
    replacement = "",
    fixed = TRUE
  )

log_debug(x = "loading feature table")
feature_table <- readr::read_delim(file = FEATURES)

log_debug(x = "loading raw files (can take long if loading multiple files)")
dda_data <- MSnbase::readMSData(
  files = FILE_POSITIVE,
  mode = "onDisk",
  msLevel. = 1
)

log_debug(x = "opening raw files objects and extracting chromatograms")
chromatograms_all <- lapply(FILE_POSITIVE, mzR::openMSfile) |>
  lapply(mzR::chromatograms) |>
  purrr::flatten()

log_debug(x = "preparing feature list ...")
set.seed(42)
df_features <- feature_table |>
  ## TODO
  tidytable::slice_sample(n = 100) |>
  prepare_features(min_intensity = 1E5)

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
  future.apply::future_lapply(
    X = list_ms_chromatograms,
    FUN = transform_ms
  )

log_debug(x = "extracting ms peaks")
list_ms_peaks <-
  future.apply::future_lapply(
    X = list_ms_chromatograms_transformed,
    FUN = extract_ms_peak
  )

log_debug(x = "comparing peaks")
list_comparison_score <-
  future.apply::future_lapply(
    X = seq_along(list_ms_peaks),
    FUN = compare_peaks
  )

log_debug(x = "selecting features with peaks")
df_features_with_peaks <-
  peaks_prelist$list_df_features_with_peaks_long |>
  tidytable::bind_rows()

log_debug(x = "There are", nrow(df_features_with_peaks), "features with peaks")

log_debug(x = "summarizing comparison scores")
comparison_scores <- list_comparison_score |>
  purrr::flatten()

log_debug(x = "There are", length(comparison_scores), "scores calculated")

log_debug(x = "joining")
df_features_with_peaks$comparison_score <-
  as.numeric(comparison_scores)

log_debug(x = "final aesthetics")
df_features_with_peaks_scored <- df_features_with_peaks |>
  tidytable::select(
    sample = id,
    peak_id,
    peak_rt_min = rt_min,
    peak_rt_apex = rt_apex,
    peak_rt_max = rt_max,
    peak_area = integral,
    feature_id,
    feature_rt = rt,
    feature_mz = mz,
    feature_area = area,
    comparison_score
  ) |>
  tidytable::distinct()

df_features_without_peaks_scored <-
  peaks_prelist$df_features_without_peaks |>
  tidytable::mutate(comparison_score = NA) |>
  tidytable::select(
    sample = id,
    peak_id,
    peak_rt_min = rt_min,
    peak_rt_apex = rt_apex,
    peak_rt_max = rt_max,
    peak_area = integral,
    feature_id,
    feature_rt = rt,
    feature_mz = mz,
    feature_area = area,
    comparison_score
  ) |>
  tidytable::distinct()

log_debug(x = "checking export directory")
check_export_dir(EXPORT_DIR)

log_debug(x = "exporting to ...")
log_debug(x = EXPORT_DIR)
readr::write_tsv(
  x = df_features_with_peaks_scored,
  file = file.path(EXPORT_DIR, switch(detector,
    "bpi" = EXPORT_FILE_BPI,
    "cad" = EXPORT_FILE_CAD,
    "pda" = EXPORT_FILE_PDA
  ))
)
## TODO FIX THIS
EXPORT_FILE_CAD_2 <- EXPORT_FILE_CAD_2[1]
readr::write_tsv(
  x = df_features_without_peaks_scored,
  file = file.path(EXPORT_DIR, switch(detector,
    "bpi" = EXPORT_FILE_BPI_2,
    "cad" = EXPORT_FILE_CAD_2,
    "pda" = EXPORT_FILE_PDA_2
  ))
)

# if (params$signal$detector$pda == TRUE) {}
# TODO check if change for PDA doing baselining on normal and not improved
# chromatograms_list_pda <-
#   preprocess_chromatograms(
#     detector = "pda",
#     list = chromatograms_all[c(FALSE, TRUE, FALSE)],
#     signal_name = "PDA.1_TotalAbsorbance_0",
#     shift = PDA_SHIFT
#   )
# peaks_prelist_pda <- preprocess_peaks(
#   detector = "pda",
#   list = chromatograms_list_pda$chromatograms_improved,
#   df_long = chromatograms_list_pda$chromatograms_improved_long
# )
# detector <- "pda"
# peaks_list_pda <- process_peaks(detector = "pda")

end <- Sys.time()

log_debug("Script finished in", format(end - start))
