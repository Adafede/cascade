#' Add a chromatogram trace to a plotly figure
#'
#' @param plot A `plotly` object.
#' @param chromato Data frame / `data.table` with `rtime` and `intensity`.
#' @param shift Retention time shift in minutes (default `0`).
#' @param normalize_time Logical; normalize time to `[0, 1]`.
#' @param name Trace display name.
#' @param color Line color hex code.
#' @param polarity Either `"pos"` or `"neg"`.
#'
#' @return Updated `plotly` object.
#' @keywords internal
#' @noRd
add_chromato_line <- function(
  plot,
  chromato,
  shift = 0,
  normalize_time = FALSE,
  name,
  color,
  polarity = "pos"
) {
  norm_df <- normalize_chromatograms_list(
    chromato,
    shift = shift,
    normalize_time = normalize_time
  )
  if (polarity == "neg") {
    norm_df <- data.table::copy(data.table::as.data.table(norm_df))
    norm_df[, intensity := -1 * intensity]
  }
  plotly::add_lines(
    plot,
    data = norm_df,
    x = ~rtime,
    y = ~intensity,
    name = name,
    line = list(width = 1, color = color)
  )
}

#' Build interactive multi-detector chromatogram comparison plot
#'
#' @param chromatograms Character vector of chromatogram keys to display.
#' @param chromatograms_list Named list of processed chromatogram tables.
#' @param normalize_time Logical; normalize retention time to `[0, 1]`.
#' @param shift_cad CAD retention time shift in minutes.
#' @param shift_pda PDA retention time shift in minutes.
#' @param type Either `"improved"` or `"baselined"`.
#'
#' @return A `plotly` object.
#' @keywords internal
#' @noRd
check_chromatograms <- function(
  chromatograms = c("bpi_pos", "cad_pos", "pda_pos"),
  chromatograms_list,
  normalize_time = FALSE,
  shift_cad = 0,
  shift_pda = 0,
  type = "improved"
) {
  valid_chroms <- c(
    "bpi_pos",
    "cad_pos",
    "pda_pos",
    "bpi_neg",
    "cad_neg",
    "pda_neg"
  )
  if (!all(chromatograms %in% valid_chroms)) {
    stop(
      "chromatograms must be in 'bpi_pos', 'bpi_neg', ",
      "'cad_pos', 'cad_neg', 'pda_pos', 'pda_neg'",
      call. = FALSE
    )
  }
  stopifnot(
    "type must be one of 'improved' or 'baselined'" = type %in%
      c("improved", "baselined")
  )

  specs <- list(
    cad_pos = list(
      key = paste0("chromatogram_cad_pos_", type),
      shift = shift_cad,
      name = "<b> CAD </b>",
      color = "#e31a1c",
      polarity = "pos"
    ),
    cad_neg = list(
      key = paste0("chromatogram_cad_neg_", type),
      shift = shift_cad,
      name = "<b> CAD (Neg)</b>",
      color = "#fb9a99",
      polarity = "neg"
    ),
    pda_pos = list(
      key = paste0("chromatogram_pda_pos_", type),
      shift = shift_pda,
      name = "<b> PDA </b>",
      color = "#33a02c",
      polarity = "pos"
    ),
    pda_neg = list(
      key = paste0("chromatogram_pda_neg_", type),
      shift = shift_pda,
      name = "<b> PDA (Neg) </b>",
      color = "#b2df8a",
      polarity = "neg"
    ),
    bpi_pos = list(
      key = paste0("chromatogram_bpi_pos_", type),
      shift = 0,
      name = "<b> MS Pos </b>",
      color = "#1f78b4",
      polarity = "pos"
    ),
    bpi_neg = list(
      key = paste0("chromatogram_bpi_neg_", type),
      shift = 0,
      name = "<b> MS Neg </b>",
      color = "#a6cee3",
      polarity = "neg"
    )
  )

  plot <- plotly::plot_ly()
  for (nm in names(specs)) {
    if (nm %in% chromatograms) {
      sp <- specs[[nm]]
      chrom_data <- chromatograms_list[[sp$key]]
      if (!is.null(chrom_data)) {
        plot <- add_chromato_line(
          plot = plot,
          chromato = chrom_data,
          shift = sp$shift,
          normalize_time = normalize_time,
          name = sp$name,
          color = sp$color,
          polarity = sp$polarity
        )
      }
    }
  }

  plotly::layout(
    plot,
    xaxis = list(
      title = if (normalize_time) {
        "<b> Normalized Time </b>"
      } else {
        "<b> Time [minutes] </b>"
      }
    ),
    yaxis = list(title = "<b> Normalized Intensity </b>")
  )
}

#' Check multi-detector chromatogram alignment
#'
#' Visualizes and aligns BPI, CAD, and PDA chromatograms across positive and/or
#' negative ionization runs after optional Fourier noise filtering, peak
#' sharpening, and baseline correction.
#'
#' @export
#'
#' @param file_negative Path to negative-mode `.mzML` file (optional).
#' @param file_positive Path to positive-mode `.mzML` file (optional).
#' @param time_min Minimum retention time in minutes. Default is `0.5`.
#' @param time_max Maximum retention time in minutes. Default is `32.5`.
#' @param cad_shift CAD retention time shift in minutes, or `"auto"` / `NULL`
#'   to estimate the optimal shift via FFT cross-correlation against the MS BPI
#'   trace. Default is `0.05`.
#' @param pda_shift PDA retention time shift in minutes, or `"auto"` / `NULL`
#'   to estimate the optimal shift via FFT cross-correlation against the MS BPI
#'   trace. Default is `0.1`.
#' @param fourier_components Fraction of Fourier components to retain. Default
#'   is `0.01`.
#' @param frequency Acquisition frequency in Hz. Default is `1`.
#' @param resample Resampling factor. Default is `1`.
#' @param chromatograms Character vector of chromatograms to plot. Default is
#'   `c("bpi_pos", "cad_pos", "pda_pos")`.
#' @param headers Named character vector mapping detector types (`"bpi"`,
#'   `"pda"`, `"cad"`) to `.mzML` chromatogram header IDs.
#' @param type Type of chromatogram to display: `"baselined"` or `"improved"`.
#'   Default is `"baselined"`.
#' @param normalize_intensity Logical; normalize intensity to `[0, 1]`. Default
#'   is `TRUE`.
#' @param normalize_time Logical; normalize time to `[0, 1]`. Default is
#'   `FALSE`.
#' @param show_example Logical; use bundled example data. Default is `FALSE`.
#' @param intensity_floor Small positive value for intensity floor. Default is
#'   `0.001`.
#' @param k2 Second-derivative weight for signal sharpening. Default is `250`.
#' @param k4 Fourth-derivative weight for signal sharpening. Default is
#'   `1250000`.
#' @param sigma Gaussian sigma for signal sharpening. Default is `0.05`.
#' @param smoothing_width Moving-average window width for derivative smoothing.
#'   Default is `8`.
#' @param baseline_method Baseline correction method passed to
#'   [baseline::baseline()]. Default is `"peakDetection"`.
#' @param improve_signal Logical; apply Fourier filtering and sharpening.
#'   Default is `TRUE`.
#'
#' @return A `plotly` htmlwidget displaying the aligned chromatograms.
#'
#' @concept EDAM:operation_3628 Chromatographic alignment
#' @concept EDAM:operation_3203 Chromatogram visualisation
#' @concept EDAM:topic_3172 Metabolomics
#'
#' @examples
#' check_chromatograms_alignment(show_example = TRUE)
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
  chromatogram_selected <- function(pattern) {
    any(grepl(pattern, chromatograms, fixed = TRUE))
  }

  process_chromatogram <- function(chromatogram_df) {
    if (improve_signal) {
      improved <- improve_signal(
        chromatogram_df,
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
      dt <- data.table::as.data.table(chromatogram_df)
      improved <- as.data.frame(
        dt[rtime >= time_min & rtime <= time_max, list(rtime, intensity)]
      )
    }

    improved_normalized <- normalize_chromatograms_list(
      improved,
      normalize_intensity = normalize_intensity,
      normalize_time = normalize_time
    )

    baselined <- normalize_chromatograms_list(
      baseline_chromatogram(improved_normalized, method = baseline_method),
      normalize_intensity = normalize_intensity,
      normalize_time = normalize_time
    )

    list(improved = improved_normalized, baselined = baselined)
  }

  chromatograms_list <- list()
  detectors <- c("bpi", "cad", "pda")

  if (!is.null(file_positive) || show_example) {
    chromatograms_positive <- if (show_example) {
      load_chromatograms(
        show_example = TRUE,
        headers = headers,
        example_polarity = "pos"
      )
    } else {
      load_chromatograms(file_positive, headers = headers)
    }

    for (det in detectors) {
      if (chromatogram_selected(det)) {
        chrom_df <- extract_chromatogram(
          chromatograms_positive,
          type = det,
          headers = headers
        )
        res <- process_chromatogram(chrom_df)
        chromatograms_list[[paste0(
          "chromatogram_",
          det,
          "_pos_improved"
        )]] <- res$improved
        chromatograms_list[[paste0(
          "chromatogram_",
          det,
          "_pos_baselined"
        )]] <- res$baselined
      }
    }
  }

  if (!is.null(file_negative) || show_example) {
    chromatograms_negative <- if (show_example) {
      load_chromatograms(
        show_example = TRUE,
        headers = headers,
        example_polarity = "neg"
      )
    } else {
      load_chromatograms(file_negative, headers = headers)
    }

    for (det in detectors) {
      if (chromatogram_selected(det)) {
        chrom_df <- extract_chromatogram(
          chromatograms_negative,
          type = det,
          headers = headers
        )
        res <- process_chromatogram(chrom_df)
        chromatograms_list[[paste0(
          "chromatogram_",
          det,
          "_neg_improved"
        )]] <- res$improved
        chromatograms_list[[paste0(
          "chromatogram_",
          det,
          "_neg_baselined"
        )]] <- res$baselined
      }
    }
  }

  ref_bpi <- chromatograms_list[[paste0("chromatogram_bpi_pos_", type)]]
  if (is.null(ref_bpi)) {
    ref_bpi <- chromatograms_list[[paste0("chromatogram_bpi_neg_", type)]]
  }
  if (is.null(cad_shift) || identical(cad_shift, "auto")) {
    tgt_cad <- chromatograms_list[[paste0("chromatogram_cad_pos_", type)]]
    if (is.null(tgt_cad)) {
      tgt_cad <- chromatograms_list[[paste0("chromatogram_cad_neg_", type)]]
    }
    cad_shift <- if (!is.null(ref_bpi) && !is.null(tgt_cad)) {
      estimate_rt_shift(ref_bpi, tgt_cad)
    } else {
      0
    }
  }
  if (is.null(pda_shift) || identical(pda_shift, "auto")) {
    tgt_pda <- chromatograms_list[[paste0("chromatogram_pda_pos_", type)]]
    if (is.null(tgt_pda)) {
      tgt_pda <- chromatograms_list[[paste0("chromatogram_pda_neg_", type)]]
    }
    pda_shift <- if (!is.null(ref_bpi) && !is.null(tgt_pda)) {
      estimate_rt_shift(ref_bpi, tgt_pda)
    } else {
      0
    }
  }

  check_chromatograms(
    chromatograms = chromatograms,
    chromatograms_list = chromatograms_list,
    normalize_time = normalize_time,
    shift_cad = cad_shift,
    shift_pda = pda_shift,
    type = type
  )
}
