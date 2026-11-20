#' Compare MS extracted ion chromatograms against detector peaks
#'
#' Extracts MS1 ion chromatograms (XICs) for all features falling within
#' detected and trapezoidally integrated chromatographic peaks (e.g., CAD, PDA,
#' or BPI), computes window-preserving interpolated Pearson correlation
#' similarity scores between each feature's XIC profile and the corresponding
#' detector peak profile, and exports both peak-informed and non-informed
#' feature tables as TSV files.
#'
#' @export
#'
#' @param file Path to `.mzML` file (optional when `show_example = TRUE`).
#' @param features Path to feature table (optional when `show_example = TRUE`).
#' @param type Chromatogram processing stage: `"original"`, `"baselined"`, or
#'   `"improved"`. Default is `"baselined"`.
#' @param detector Detector type: `"cad"`, `"bpi"`, or `"pda"`. Default is
#'   `"cad"`.
#' @param headers Named character vector mapping detector types to `.mzML`
#'   chromatogram header IDs.
#' @param export_dir Directory where the output `.tsv` files are written.
#'   Default is `"data/interim/peaks"`.
#' @param show_example Logical; use bundled example data. Default is `FALSE`.
#' @param fourier_components Fraction of Fourier components to retain. Default
#'   is `0.01`.
#' @param frequency Acquisition frequency in Hz. Default is `1`.
#' @param min_area Minimum relative area fraction for peak filtering. Default is
#'   `0.005`.
#' @param min_intensity Minimum MS intensity for feature filtering. Default is
#'   `1e4`.
#' @param resample Resampling factor. Default is `1`.
#' @param shift Retention time shift in minutes. Default is `0.05`.
#' @param time_min Minimum retention time in minutes. Default is `0.5`.
#' @param time_max Maximum retention time in minutes. Default is `32.5`.
#'
#' @return Invisibly returns a named list with two `data.table` elements:
#'   `features_informed` and `features_not_informed`.
#'
#' @concept EDAM:operation_3215 Peak detection
#' @concept EDAM:operation_3628 Chromatographic alignment
#' @concept EDAM:operation_3799 Quantification
#' @concept EDAM:data_2536 Mass spectrometry data
#' @concept EDAM:topic_3172 Metabolomics
#'
#' @examples
#' \donttest{
#' tmp <- tempdir()
#' process_compare_peaks(show_example = TRUE, export_dir = tmp)
#' }
process_compare_peaks <- function(
  file = NULL,
  features = NULL,
  type = "baselined",
  detector = "cad",
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  export_dir = "data/interim/peaks",
  show_example = FALSE,
  fourier_components = 0.01,
  frequency = 1,
  min_area = 0.005,
  min_intensity = 1E4,
  resample = 1,
  shift = 0.05,
  time_min = 0.5,
  time_max = 32.5
) {
  detector <- match.arg(detector, c("cad", "bpi", "pda"))
  type <- match.arg(type, c("baselined", "improved", "original"))

  message("loading MS data")
  ms_data <- load_ms_data(file, show_example = show_example)

  chrom_sub <- load_detector_chromatogram(
    file = file,
    detector = detector,
    headers = headers,
    show_example = show_example
  )

  message("loading name")
  name <- load_name(file, show_example = show_example)

  message("loading feature table")
  feature_table <- load_features(features, show_example = show_example)

  message("preparing features")
  df_features <- prepare_features(
    feature_table,
    min_intensity = min_intensity,
    name = name
  )

  message("preprocessing chromatograms")
  chromatograms_list <- preprocess_chromatograms(
    detector = detector,
    name = name,
    list = chrom_sub,
    shift = 0,
    fourier_components = fourier_components,
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample
  )

  stage <- select_chromatogram_stage(chromatograms_list, type)
  df_long_sel <- stage$df_long
  df_xy_sel <- stage$df_xy

  message("preprocessing peaks")
  peaks_prelist <- preprocess_peaks(
    detector = detector,
    df_features = df_features,
    df_long = df_long_sel,
    df_xy = df_xy_sel,
    min_area = min_area,
    shift = shift,
    name = name
  )

  message("processing ", detector, " peaks")
  message("extracting ms chromatograms (longest step)")
  n_groups <- length(peaks_prelist$list_df_features_with_peaks_long)
  list_ms_chromatograms <- extract_ms_progress(
    xs = seq_len(n_groups),
    ms_data = ms_data,
    rts = peaks_prelist$list_rtr,
    mzs = peaks_prelist$list_mzr
  )

  message("extracting ms peaks")
  list_ms_peaks <- lapply(list_ms_chromatograms, extract_ms_peak)

  message("comparing peaks")
  list_comparison_score <- lapply(
    seq_along(list_ms_peaks),
    compare_peaks,
    list_ms_peaks = list_ms_peaks,
    peaks_prelist = peaks_prelist
  )

  message("summarizing comparison scores")
  n_scores <- sum(lengths(list_comparison_score))
  message("there are ", n_scores, " scores calculated")

  message("selecting features with peaks")
  list_df_features_with_scores <- Map(
    function(df, score) {
      df[, comparison_score := unlist(score, use.names = FALSE)]
      df
    },
    peaks_prelist$list_df_features_with_peaks_long,
    list_comparison_score
  )

  if (length(list_df_features_with_scores) == 0L) {
    df_features_with_peaks <- peaks_prelist$df_features_without_peaks[0L, ]
    df_features_with_peaks[, comparison_score := numeric(0)]
  } else {
    df_features_with_peaks <- data.table::rbindlist(
      list_df_features_with_scores,
      use.names = TRUE,
      fill = TRUE
    )
  }

  message("final aesthetics")
  format_scored_features <- function(dt) {
    unique(dt[, list(
      sample = id,
      peak_id = peak_id,
      peak_rt_min = rt_min,
      peak_rt_apex = rt_apex,
      peak_rt_max = rt_max,
      peak_area = integral,
      feature_id = feature_id,
      feature_rt = rt,
      feature_mz = mz,
      feature_area = area,
      comparison_score = comparison_score
    )])
  }

  df_features_with_peaks_scored <- format_scored_features(
    df_features_with_peaks
  )

  without_dt <- peaks_prelist$df_features_without_peaks
  without_dt[, comparison_score := NA_real_]
  df_feat_without_peaks_scored <- format_scored_features(without_dt)

  message("checking export directory")
  check_export_dir(export_dir)

  message("exporting")
  base_stem <- sub("\\.[^.]+$", "", name)
  data.table::fwrite(
    df_features_with_peaks_scored,
    file = file.path(
      export_dir,
      paste0(base_stem, "_featuresInformed_", detector, ".tsv")
    ),
    sep = "\t"
  )

  data.table::fwrite(
    df_feat_without_peaks_scored,
    file = file.path(
      export_dir,
      paste0(base_stem, "_featuresNotInformed_", detector, ".tsv")
    ),
    sep = "\t"
  )

  invisible(list(
    features_informed = df_features_with_peaks_scored,
    features_not_informed = df_feat_without_peaks_scored
  ))
}
