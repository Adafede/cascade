#' Plot chromatogram with detected peak apexes and integration bounds
#'
#' @param df1 Chromatogram `data.table` with `rtime` and `intensity`.
#' @param df2 Peaks `data.table` with `rt_apex`, `peak_max`, `rt_min`, `rt_max`.
#' @param fun Approximating function mapping retention time to normalized signal
#'   intensity.
#'
#' @return A `plotly` object.
#' @keywords internal
#' @noRd
plot_peak_detection <- function(df1, df2, fun) {
  plotly::plot_ly(df1) |>
    plotly::add_trace(
      data = df1,
      x = ~rtime,
      y = ~intensity,
      type = "scatter",
      mode = "lines",
      name = "signal",
      line = list(color = "#1f78b4", width = 1)
    ) |>
    plotly::add_trace(
      data = df2,
      x = ~rt_apex,
      y = ~peak_max,
      yaxis = "y2",
      type = "scatter",
      mode = "markers",
      marker = list(color = "#ff7f00", symbol = "star"),
      name = "detected maximum",
      line = list(color = "#1f78b4", width = 0)
    ) |>
    plotly::add_trace(
      data = df2,
      x = ~rt_min,
      y = ~ fun(rt_min),
      yaxis = "y2",
      type = "scatter",
      mode = "markers",
      marker = list(color = "#33a02c", symbol = "triangle-right"),
      name = "detected minimum (start)",
      line = list(color = "#1f78b4", width = 0)
    ) |>
    plotly::add_trace(
      data = df2,
      x = ~rt_max,
      y = ~ fun(rt_max),
      yaxis = "y2",
      type = "scatter",
      mode = "markers",
      marker = list(color = "#33a02c", symbol = "triangle-left"),
      name = "detected minimum (end)",
      line = list(color = "#1f78b4", width = 0)
    ) |>
    plotly::layout(
      yaxis = list(
        title = "Normalized Intensity",
        zeroline = TRUE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis2 = list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        overlaying = "y",
        side = "right"
      ),
      xaxis = list(
        title = "Time [min]",
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      showlegend = FALSE
    )
}

#' Check chromatographic peak detection and integration bounds
#'
#' Detects, fits (via Exponential-Gaussian Hybrid or Gaussian nonlinear least
#' squares), and trapezoidally integrates peaks in a detector chromatogram (CAD,
#' BPI, or PDA), links them to MS features, and returns an interactive `plotly`
#' visualization showing the signal, detected peak maxima, and start/end
#' integration bounds.
#'
#' @export
#'
#' @param file Path to `.mzML` file (optional when `show_example = TRUE`).
#' @param features Path to feature table (optional when `show_example = TRUE`).
#' @param detector Detector type: `"cad"`, `"bpi"`, or `"pda"`. Default is
#'   `"cad"`.
#' @param chromatogram Chromatogram processing stage: `"original"`,
#'   `"improved"`, or `"baselined"`. Default is `"baselined"`.
#' @param headers Named character vector mapping detector types to `.mzML`
#'   chromatogram header IDs.
#' @param min_area Minimum relative area fraction for peak filtering. Default is
#'   `0.005`.
#' @param min_intensity Minimum MS intensity for feature filtering. Default is
#'   `1e4`.
#' @param shift Retention time shift in minutes. Default is `0.05`.
#' @param show_example Logical; use bundled example data. Default is `FALSE`.
#' @param fourier_components Fraction of Fourier components to retain. Default
#'   is `0.01`.
#' @param time_min Minimum retention time in minutes. Default is `0.5`.
#' @param time_max Maximum retention time in minutes. Default is `32.5`.
#' @param frequency Acquisition frequency in Hz. Default is `1`.
#' @param resample Resampling factor. Default is `1`.
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
#' @param sd_max Maximum standard deviation for peak filtering. Default is `50`.
#' @param max_iter Maximum iterations for peak fitting. Default is `1000`.
#' @param noise_threshold Noise threshold for peak detection. Default is
#'   `0.001`.
#' @param fit Peak fitting model: `"egh"`, `"gaussian"`, or `"raw"`. Default is
#'   `"egh"`.
#' @param intensity_threshold Minimum normalized intensity threshold for peak
#'   filtering. Default is `0.1`.
#' @param improve_signal Logical; apply Fourier filtering and sharpening.
#'   Default is `TRUE`.
#'
#' @return A `plotly` htmlwidget showing the chromatogram and detected peak
#'   boundaries.
#'
#' @concept EDAM:operation_3215 Peak detection
#' @concept EDAM:operation_3203 Chromatogram visualisation
#' @concept EDAM:topic_3172 Metabolomics
#'
#' @examples
#' check_peaks_integration(show_example = TRUE)
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
  intensity_threshold = 0.1,
  improve_signal = TRUE
) {
  detector <- match.arg(detector, c("cad", "bpi", "pda"))
  chromatogram <- match.arg(
    chromatogram,
    c("baselined", "improved", "original")
  )

  message("loading chromatograms")
  chromatograms_all <- load_chromatograms(
    file,
    show_example = show_example,
    headers = headers
  )
  if (
    !(detector %in% names(chromatograms_all)) ||
      is.null(chromatograms_all[[detector]])
  ) {
    stop(
      "Detector '",
      detector,
      "' chromatogram not found in file.",
      call. = FALSE
    )
  }

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

  message("Preprocessing chromatograms")
  chrom_sub <- chromatograms_all[detector]
  chromatograms_list <- preprocess_chromatograms(
    detector = detector,
    name = name,
    list = chrom_sub,
    shift = 0,
    fourier_components = fourier_components,
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample,
    intensity_floor = intensity_floor,
    k2 = k2,
    k4 = k4,
    sigma = sigma,
    smoothing_width = smoothing_width,
    baseline_method = baseline_method,
    improve_signal = improve_signal
  )

  df_long_sel <- data.table::copy(data.table::as.data.table(switch(
    chromatogram,
    "original" = chromatograms_list$chromatograms_original_long,
    "improved" = chromatograms_list$chromatograms_improved_long,
    "baselined" = chromatograms_list$chromatograms_baselined_long
  )))
  max_int <- max(df_long_sel$intensity, na.rm = TRUE)
  if (is.finite(max_int) && max_int > 0) {
    df_long_sel[, intensity := intensity / max_int]
  }

  df_xy_sel <- switch(
    chromatogram,
    "original" = chromatograms_list$chromatograms_original[[1L]],
    "improved" = chromatograms_list$chromatograms_improved[[1L]],
    "baselined" = chromatograms_list$chromatograms_baselined[[1L]]
  )

  peaks <- preprocess_peaks(
    df_features = df_features,
    df_long = df_long_sel,
    df_xy = df_xy_sel,
    min_area = min_area,
    shift = shift,
    name = name,
    sd_max = sd_max,
    max_iter = max_iter,
    noise_threshold = noise_threshold,
    fit = fit,
    intensity_threshold = intensity_threshold
  )

  chromatogram_normalized <- if (chromatogram == "original") {
    orig_dt <- data.table::copy(data.table::as.data.table(
      chromatograms_list$chromatograms_original_long
    ))
    orig_dt <- orig_dt[seq_len(nrow(orig_dt)) %% 10L == 1L]
    max_orig <- max(orig_dt$intensity, na.rm = TRUE)
    if (is.finite(max_orig) && max_orig > 0) {
      orig_dt[, intensity := intensity / max_orig]
    }
    orig_dt
  } else {
    df_long_sel
  }

  peaks_normalized <- data.table::rbindlist(
    peaks$list_df_features_with_peaks_long,
    use.names = TRUE,
    fill = TRUE
  )
  if (nrow(peaks_normalized) > 0L) {
    max_imax <- max(peaks_normalized$intensity_max, na.rm = TRUE)
    max_pmax <- max(peaks_normalized$peak_max, na.rm = TRUE)
    if (is.finite(max_imax) && max_imax > 0) {
      peaks_normalized[, intensity := intensity_max / max_imax]
    } else {
      peaks_normalized[, intensity := intensity_max]
    }
    if (is.finite(max_pmax) && max_pmax > 0) {
      peaks_normalized[, peak_max := peak_max / max_pmax]
    }
  }

  approx_f <- stats::approxfun(
    x = chromatogram_normalized$rtime,
    y = chromatogram_normalized$intensity
  )

  plot_peak_detection(
    df1 = chromatogram_normalized,
    df2 = peaks_normalized,
    fun = approx_f
  )
}
