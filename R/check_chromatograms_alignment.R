#' Check chromatograms alignment
#'
#' @export
#'
#' @include baseline_chromatogram.R
#' @include check_chromatograms.R
#' @include extract_chromatogram.R
#' @include improve_signal.R
#' @include load_chromatograms.R
#' @include normalize_chromatograms_list.R
#'
#' @param file_negative Negative file path
#' @param file_positive Positive file path
#' @param time_min Minimum time in minutes. Default is 0.5.
#' @param time_max Maximum time in minutes. Default is 32.5.
#' @param cad_shift CAD time shift in minutes. Default is 0.05.
#' @param pda_shift PDA time shift in minutes. Default is 0.1.
#' @param fourier_components Fraction of Fourier components to keep. Default is
#'   0.01.
#' @param frequency Acquisition frequency in Hz. Default is 1.
#' @param resample Resampling factor. Default is 1.
#' @param chromatograms Chromatograms to plot. Default is c("bpi_pos",
#'   "cad_pos", "pda_pos").
#' @param headers Named vector mapping detector types to header names in the
#'   mzML file.
#' @param type Type of chromatogram to display. Either "baselined" or
#'   "improved". Default is "baselined".
#' @param normalize_intensity Normalize intensity? Default is TRUE.
#' @param normalize_time Normalize time? Default is FALSE.
#' @param show_example Show example data? Default is FALSE.
#' @param intensity_floor Small positive value for intensity floor. Default is
#'   0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#' @param baseline_method Method for baseline correction. Default is
#'   "peakDetection".
#' @param improve_signal Logical. Whether to apply signal improvement (Fourier
#'   filtering and sharpening). Default is TRUE.
#'
#' @return A plot with (non-)aligned chromatograms
#'
#' @examples
#' \dontrun{
#' check_chromatograms_alignment(show_example = TRUE)
#' }
check_chromatograms_alignment <- function(
  file_negative = NULL,
  file_positive = NULL,
  time_min = 0.5,
  time_max = 32.5,
  cad_shift = 0.05,
  pda_shift = 0.1,
  fourier_components = 0.01,
  frequency = 1,
  resample = 1,
  chromatograms = c("bpi_pos", "cad_pos", "pda_pos"),
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  type = "baselined",
  normalize_intensity = TRUE,
  normalize_time = FALSE,
  show_example = FALSE,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8,
  baseline_method = "peakDetection",
  improve_signal = TRUE
) {
  ## Helper function to process a single chromatogram
  process_chromatogram <- function(
    chromatogram_df,
    improve_signal,
    time_min,
    time_max,
    fourier_components,
    frequency,
    resample,
    intensity_floor,
    k2,
    k4,
    sigma,
    smoothing_width,
    normalize_intensity,
    normalize_time,
    baseline_method
  ) {
    if (improve_signal) {
      improved <- chromatogram_df |>
        improve_signal(
          time_min = time_min,
          time_max = time_max,
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample,
          intensity_floor = intensity_floor,
          k2 = k2,
          k4 = k4,
          sigma = sigma,
          smoothing_width = smoothing_width
        )
    } else {
      ## Skip signal improvement - just filter by time range
      improved <- chromatogram_df |>
        tidytable::filter(rtime >= time_min & rtime <= time_max) |>
        tidytable::select(rtime, intensity) |>
        data.frame()
    }

    improved_normalized <- improved |>
      normalize_chromatograms_list(
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time
      )

    baselined <- improved_normalized |>
      baseline_chromatogram(method = baseline_method) |>
      normalize_chromatograms_list(
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time
      )

    list(improved = improved_normalized, baselined = baselined)
  }

  chromatograms_list <- list()
  if (!is.null(file_positive) || show_example) {
    if (show_example) {
      chromatograms_positive <- load_chromatograms(
        show_example = show_example,
        headers = headers,
        example_polarity = "pos"
      )
    } else {
      chromatograms_positive <- file_positive |>
        load_chromatograms(headers = headers)
    }

    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "bpi") |>
        any()
    ) {
      chromatogram_bpi_pos <- chromatograms_positive |>
        extract_chromatogram(type = "bpi", headers = headers)
      result <- process_chromatogram(
        chromatogram_df = chromatogram_bpi_pos,
        improve_signal = improve_signal,
        time_min = time_min,
        time_max = time_max,
        fourier_components = fourier_components,
        frequency = frequency,
        resample = resample,
        intensity_floor = intensity_floor,
        k2 = k2,
        k4 = k4,
        sigma = sigma,
        smoothing_width = smoothing_width,
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time,
        baseline_method = baseline_method
      )
      chromatograms_list$chromatogram_bpi_pos_improved <- result$improved
      chromatograms_list$chromatogram_bpi_pos_baselined <- result$baselined
    }
    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "cad") |>
        any()
    ) {
      chromatogram_cad_pos <- chromatograms_positive |>
        extract_chromatogram(type = "cad", headers = headers)
      result <- process_chromatogram(
        chromatogram_df = chromatogram_cad_pos,
        improve_signal = improve_signal,
        time_min = time_min,
        time_max = time_max,
        fourier_components = fourier_components,
        frequency = frequency,
        resample = resample,
        intensity_floor = intensity_floor,
        k2 = k2,
        k4 = k4,
        sigma = sigma,
        smoothing_width = smoothing_width,
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time,
        baseline_method = baseline_method
      )
      chromatograms_list$chromatogram_cad_pos_improved <- result$improved
      chromatograms_list$chromatogram_cad_pos_baselined <- result$baselined
    }
    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "pda") |>
        any()
    ) {
      chromatogram_pda_pos <- chromatograms_positive |>
        extract_chromatogram(type = "pda", headers = headers)
      result <- process_chromatogram(
        chromatogram_df = chromatogram_pda_pos,
        improve_signal = improve_signal,
        time_min = time_min,
        time_max = time_max,
        fourier_components = fourier_components,
        frequency = frequency,
        resample = resample,
        intensity_floor = intensity_floor,
        k2 = k2,
        k4 = k4,
        sigma = sigma,
        smoothing_width = smoothing_width,
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time,
        baseline_method = baseline_method
      )
      chromatograms_list$chromatogram_pda_pos_improved <- result$improved
      chromatograms_list$chromatogram_pda_pos_baselined <- result$baselined
    }
  }

  if (!is.null(file_negative) || show_example) {
    if (show_example) {
      chromatograms_negative <- load_chromatograms(
        show_example = show_example,
        headers = headers,
        example_polarity = "neg"
      )
    } else {
      chromatograms_negative <- file_negative |>
        load_chromatograms(headers = headers)
    }

    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "bpi") |>
        any()
    ) {
      chromatogram_bpi_neg <- chromatograms_negative |>
        extract_chromatogram("bpi", headers = headers)
      result <- process_chromatogram(
        chromatogram_df = chromatogram_bpi_neg,
        improve_signal = improve_signal,
        time_min = time_min,
        time_max = time_max,
        fourier_components = fourier_components,
        frequency = frequency,
        resample = resample,
        intensity_floor = intensity_floor,
        k2 = k2,
        k4 = k4,
        sigma = sigma,
        smoothing_width = smoothing_width,
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time,
        baseline_method = baseline_method
      )
      chromatograms_list$chromatogram_bpi_neg_improved <- result$improved
      chromatograms_list$chromatogram_bpi_neg_baselined <- result$baselined
    }
    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "cad") |>
        any()
    ) {
      chromatogram_cad_neg <- chromatograms_negative |>
        extract_chromatogram("cad", headers = headers)
      result <- process_chromatogram(
        chromatogram_df = chromatogram_cad_neg,
        improve_signal = improve_signal,
        time_min = time_min,
        time_max = time_max,
        fourier_components = fourier_components,
        frequency = frequency,
        resample = resample,
        intensity_floor = intensity_floor,
        k2 = k2,
        k4 = k4,
        sigma = sigma,
        smoothing_width = smoothing_width,
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time,
        baseline_method = baseline_method
      )
      chromatograms_list$chromatogram_cad_neg_improved <- result$improved
      chromatograms_list$chromatogram_cad_neg_baselined <- result$baselined
    }
    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "pda") |>
        any()
    ) {
      chromatogram_pda_neg <- chromatograms_negative |>
        extract_chromatogram("pda", headers = headers)
      result <- process_chromatogram(
        chromatogram_df = chromatogram_pda_neg,
        improve_signal = improve_signal,
        time_min = time_min,
        time_max = time_max,
        fourier_components = fourier_components,
        frequency = frequency,
        resample = resample,
        intensity_floor = intensity_floor,
        k2 = k2,
        k4 = k4,
        sigma = sigma,
        smoothing_width = smoothing_width,
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time,
        baseline_method = baseline_method
      )
      chromatograms_list$chromatogram_pda_neg_improved <- result$improved
      chromatograms_list$chromatogram_pda_neg_baselined <- result$baselined
    }
  }

  check_chromatograms(
    chromatograms = chromatograms,
    chromatograms_list = chromatograms_list,
    shift_cad = cad_shift,
    shift_pda = pda_shift,
    type = type
  )
}
