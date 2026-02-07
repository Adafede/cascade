#' Check chromatograms alignment
#'
#' @export
#'
#' @include load_chromatograms.R
#' @include load_features.R
#' @include load_name.R
#' @include plot_peak_detection.R
#' @include prepare_features.R
#' @include preprocess_chromatograms.R
#' @include preprocess_peaks.R
#'
#' @param file File path
#' @param features Features path
#' @param detector Detector
#' @param chromatogram Chromatogram
#' @param headers Headers
#' @param min_area Minimum area
#' @param min_intensity Minimum intensity
#' @param shift shift
#' @param show_example Show example? Default to FALSE
#' @param fourier_components Fourier components
#' @param time_min Time min
#' @param time_max Time max
#' @param frequency Frequency
#' @param resample Resample
#' @param intensity_offset Offset to add to intensity values to handle negative
#'   intensities. Default is 100.
#' @param intensity_floor Small value subtracted from minimum intensity.
#'   Default is 0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#' @param baseline_method Method for baseline correction. Default is
#'   "peakDetection".
#' @param sd_max Maximum standard deviation for peak filtering. Default is 50.
#' @param max_iter Maximum iterations for peak fitting. Default is 1000.
#' @param noise_threshold Noise threshold for peak detection. Default is 0.001.
#' @param fit Peak fitting method. One of "egh", "gaussian", or "raw". Default
#'   is "egh".
#' @param intensity_threshold Minimum normalized intensity threshold for
#'   filtering. Default is 0.1.
#'
#' @return A plot with (non-)aligned chromatograms
#'
#' @examples
#' \dontrun{
#' check_peaks_integration(show_example = TRUE)
#' }
check_peaks_integration <- function(
  file = NULL,
  features = NULL,
  detector = "cad",
  chromatogram = "baselined",
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  min_area = 0.005,
  min_intensity = 1E4,
  shift = 0.05,
  show_example = FALSE,
  fourier_components = 0.01,
  time_min = 0.5,
  time_max = 32.5,
  frequency = 1,
  resample = 1,
  intensity_offset = 100,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8,
  baseline_method = "peakDetection",
  sd_max = 50,
  max_iter = 1000,
  noise_threshold = 0.001,
  fit = "egh",
  intensity_threshold = 0.1
) {
  message("loading chromatograms")
  chromatograms_all <- file |>
    load_chromatograms(show_example = show_example, headers = headers)

  message("loading name")
  name <- file |>
    load_name(show_example = show_example)

  message("loading feature table")
  feature_table <- features |>
    load_features(show_example = show_example)

  message("preparing features")
  df_features <- feature_table |>
    prepare_features(min_intensity = min_intensity, name = name)

  message("Preprocessing chromatograms")
  switch <- switch(
    detector,
    "bpi" = headers["bpi"],
    "cad" = headers["cad"],
    "pda" = headers["pda"]
  )
  list <- chromatograms_all[switch |> names()]
  chromatograms_list <- preprocess_chromatograms(
    detector = detector,
    name = name,
    list = list,
    # signal_name = signal_name,
    shift = shift,
    fourier_components = fourier_components,
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample,
    intensity_offset = intensity_offset,
    intensity_floor = intensity_floor,
    k2 = k2,
    k4 = k4,
    sigma = sigma,
    smoothing_width = smoothing_width,
    baseline_method = baseline_method
  )

  peaks <-
    preprocess_peaks(
      df_features = df_features,
      df_long = switch(
        chromatogram,
        "original" = chromatograms_list$chromatograms_original_long,
        "improved" = chromatograms_list$chromatograms_improved_long,
        "baselined" = chromatograms_list$chromatograms_baselined_long
      ) |>
        tidytable::mutate(intensity = intensity / max(intensity)),
      df_xy = switch(
        chromatogram,
        "original" = chromatograms_list$chromatograms_original[[1]],
        "improved" = chromatograms_list$chromatograms_improved[[1]],
        "baselined" = chromatograms_list$chromatograms_baselined[[1]]
      ),
      min_area = min_area,
      shift = shift,
      name = name,
      sd_max = sd_max,
      max_iter = max_iter,
      noise_threshold = noise_threshold,
      fit = fit,
      intensity_threshold = intensity_threshold
    )

  chromatogram_normalized <- switch(
    chromatogram,
    "original" = chromatograms_list$chromatograms_original_long |>
      tidytable::bind_rows() |>
      tidytable::filter(row_number() %% 10 == 1) |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    "improved" = chromatograms_list$chromatograms_improved_long |>
      tidytable::bind_rows() |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    "baselined" = chromatograms_list$chromatograms_baselined_long |>
      tidytable::bind_rows() |>
      tidytable::mutate(intensity = intensity / max(intensity))
  )

  peaks_normalized <- peaks$list_df_features_with_peaks_long |>
    tidytable::bind_rows() |>
    tidytable::mutate(
      intensity = intensity_max / max(intensity_max),
      peak_max = peak_max / max(peak_max)
    )

  approx_f <- stats::approxfun(
    x = chromatogram_normalized |>
      tidytable::pull(rtime),
    y = chromatogram_normalized |>
      tidytable::pull(intensity)
  )

  chromatogram_normalized |>
    plot_peak_detection(df2 = peaks_normalized, fun = approx_f)
}
