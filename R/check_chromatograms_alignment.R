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
#' @param time_min Minimum time
#' @param time_max Maximum time
#' @param cad_shift CAD shift
#' @param pda_shift PDA shift
#' @param fourier_components Fourier components
#' @param frequency Frequency
#' @param resample Resample
#' @param chromatograms Chromatograms to plot
#' @param headers Headers
#' @param type Type. "baselined" or "improved"
#' @param normalize_intensity Normalize intensity? Default to TRUE
#' @param normalize_time Normalize time? Default to FALSE
#' @param show_example Show example? Default to FALSE
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
  show_example = FALSE
) {
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
      chromatograms_list$chromatogram_bpi_pos_improved <- chromatogram_bpi_pos |>
        improve_signal(
          time_min = time_min,
          time_max = time_max,
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample
        ) |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
      chromatograms_list$chromatogram_bpi_pos_baselined <- chromatograms_list$chromatogram_bpi_pos_improved |>
        baseline_chromatogram() |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
    }
    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "cad") |>
        any()
    ) {
      chromatogram_cad_pos <- chromatograms_positive |>
        extract_chromatogram(type = "cad", headers = headers)
      chromatograms_list$chromatogram_cad_pos_improved <- chromatogram_cad_pos |>
        improve_signal(
          time_min = time_min,
          time_max = time_max,
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample
        ) |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
      chromatograms_list$chromatogram_cad_pos_baselined <- chromatograms_list$chromatogram_cad_pos_improved |>
        baseline_chromatogram() |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
    }
    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "pda") |>
        any()
    ) {
      chromatogram_pda_pos <- chromatograms_positive |>
        extract_chromatogram(type = "pda", headers = headers)
      chromatograms_list$chromatogram_pda_pos_improved <- chromatogram_pda_pos |>
        improve_signal(
          time_min = time_min,
          time_max = time_max,
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample
        ) |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
      chromatograms_list$chromatogram_pda_pos_baselined <- chromatograms_list$chromatogram_pda_pos_improved |>
        baseline_chromatogram() |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
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
      chromatograms_list$chromatogram_bpi_neg_improved <- chromatogram_bpi_neg |>
        improve_signal(
          time_min = time_min,
          time_max = time_max,
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample
        ) |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
      chromatograms_list$chromatogram_bpi_neg_baselined <- chromatograms_list$chromatogram_bpi_neg_improved |>
        baseline_chromatogram() |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
    }
    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "cad") |>
        any()
    ) {
      chromatogram_cad_neg <- chromatograms_negative |>
        extract_chromatogram("cad", headers = headers)
      chromatograms_list$chromatogram_cad_neg_improved <- chromatogram_cad_neg |>
        improve_signal(
          time_min = time_min,
          time_max = time_max,
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample
        ) |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
      chromatograms_list$chromatogram_cad_neg_baselined <- chromatograms_list$chromatogram_cad_neg_improved |>
        baseline_chromatogram() |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
    }
    if (
      chromatograms |>
        stringi::stri_detect_fixed(pattern = "pda") |>
        any()
    ) {
      chromatogram_pda_neg <- chromatograms_negative |>
        extract_chromatogram("pda", headers = headers)
      chromatograms_list$chromatogram_pda_neg_improved <- chromatogram_pda_neg |>
        improve_signal(
          time_min = time_min,
          time_max = time_max,
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample
        ) |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
      chromatograms_list$chromatogram_pda_neg_baselined <- chromatograms_list$chromatogram_pda_neg_improved |>
        baseline_chromatogram() |>
        normalize_chromatograms_list(
          normalize_intensity = normalize_intensity,
          normalize_time = normalize_time
        )
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
