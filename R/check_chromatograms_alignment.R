#' Check chromatograms alignment
#'
#' @export
#'
#' @include load_chromatograms.R
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
#' @param type Type. "baselined" or "improved"
#' @param normalize_intensity Normalize intensity? Default to TRUE
#' @param normalize_time Normalize time? Default to FALSE
#' @param show_example Show example? Default to FALSE
#'
#' @return A plot with (non-)aligned chromatograms
#'
#' @examples
#' \dontrun{
#' cascade:::check_chromatograms_alignment()
#' }
check_chromatograms_alignment <- function(file_negative = NULL,
                                          file_positive = NULL,
                                          time_min = 0.7,
                                          time_max = 35.2,
                                          cad_shift = 0.05,
                                          pda_shift = 0.1,
                                          fourier_components = 0.01,
                                          frequency = 2,
                                          resample = 1,
                                          chromatograms = c("bpi_pos", "cad_pos", "pda_pos"),
                                          type = "baselined",
                                          normalize_intensity = TRUE,
                                          normalize_time = FALSE,
                                          show_example = FALSE) {
  if (!is.null(file_positive) | show_example) {
    if (show_example) {
      chromatograms_positive <- load_chromatograms(show_example = show_example, example_polarity = "pos")
    } else {
      chromatograms_positive <- file_positive |>
        load_chromatograms()
    }

    chromatogram_bpi_pos <- chromatograms_positive |>
      extract_chromatogram("bpi")
    chromatogram_cad_pos <- chromatograms_positive |>
      extract_chromatogram("cad")
    chromatogram_pda_pos <- chromatograms_positive |>
      extract_chromatogram("pda")

    chromatogram_bpi_pos_improved <<- chromatogram_bpi_pos |>
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
    chromatogram_cad_pos_improved <<- chromatogram_cad_pos |>
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
    chromatogram_pda_pos_improved <<- chromatogram_pda_pos |>
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

    chromatogram_bpi_pos_baselined <<- chromatogram_bpi_pos_improved |>
      baseline_chromatogram() |>
      normalize_chromatograms_list(
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time
      )
    chromatogram_cad_pos_baselined <<- chromatogram_cad_pos_improved |>
      baseline_chromatogram() |>
      normalize_chromatograms_list(
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time
      )
    chromatogram_pda_pos_baselined <<- chromatogram_pda_pos_improved |>
      baseline_chromatogram() |>
      normalize_chromatograms_list(
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time
      )
  }

  if (!is.null(file_negative) | show_example) {
    if (show_example) {
      chromatograms_negative <- load_chromatograms(show_example = show_example, example_polarity = "neg")
    } else {
      chromatograms_negative <- file_negative |>
        load_chromatograms()
    }

    chromatogram_bpi_neg <- chromatograms_negative |>
      extract_chromatogram("bpi")
    chromatogram_cad_neg <- chromatograms_negative |>
      extract_chromatogram("cad")
    chromatogram_pda_neg <- chromatograms_negative |>
      extract_chromatogram("pda")

    chromatogram_bpi_neg_improved <<- chromatogram_bpi_neg |>
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
    chromatogram_cad_neg_improved <<- chromatogram_cad_neg |>
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
    chromatogram_pda_neg_improved <<- chromatogram_pda_neg |>
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

    chromatogram_bpi_neg_baselined <<- chromatogram_bpi_neg_improved |>
      baseline_chromatogram() |>
      normalize_chromatograms_list(
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time
      )
    chromatogram_cad_neg_baselined <<- chromatogram_cad_neg_improved |>
      baseline_chromatogram() |>
      normalize_chromatograms_list(
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time
      )
    chromatogram_pda_neg_baselined <<- chromatogram_pda_neg_improved |>
      baseline_chromatogram() |>
      normalize_chromatograms_list(
        normalize_intensity = normalize_intensity,
        normalize_time = normalize_time
      )
  }

  check_chromatograms(
    shift_cad = cad_shift,
    shift_pda = pda_shift,
    type = type,
    chromatograms = chromatograms
  )
}
