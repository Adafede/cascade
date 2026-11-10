#' Centered running mean matching caTools::runmean(x, k, align = "center")
#'
#' @param x Numeric vector.
#' @param k Integer(1) window width.
#'
#' @return Numeric vector of centered running means of length `length(x)`.
#' @keywords internal
#' @noRd
runmean_center <- function(x, k) {
  x <- as.vector(x)
  n <- length(x)
  k <- as.integer(k)
  if (k <= 1L || n == 0L) {
    return(x)
  }
  if (k > n) {
    k <- n
  }
  k2 <- k %/% 2L
  k1 <- k - k2 - 1L
  y <- c(
    rep(0, k1),
    cumsum(c(sum(x[seq_len(k)]), diff(x, k))) / k,
    rep(0, k2)
  )
  if (k == 2L) {
    y[2:n] <- y[1:(n - 1L)]
    y[1L] <- mean(x[1L], na.rm = TRUE)
    return(y)
  }
  if (k1 > 0L) {
    for (i in seq_len(k1)) {
      y[i] <- mean(x[seq_len(i + k2)], na.rm = TRUE)
    }
  }
  for (i in (n - k2 + 1L):n) {
    y[i] <- mean(x[(i - k1):n], na.rm = TRUE)
  }
  y
}

#' First derivative
#'
#' @param x Numeric vector of abscissa values.
#' @param y Numeric vector of ordinate values.
#'
#' @return Numeric vector of first finite differences `diff(y) / diff(x)`.
#' @keywords internal
#' @noRd
deriv <- function(x, y) {
  diff(y) / diff(x)
}

#' Midpoints of a numeric vector
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of midpoints between consecutive elements of `x`.
#' @keywords internal
#' @noRd
middle_pts <- function(x) {
  x[-1L] - diff(x) / 2
}

#' Second derivative
#'
#' @param x Numeric vector of abscissa values.
#' @param y Numeric vector of ordinate values.
#'
#' @return Numeric vector of second finite differences.
#' @keywords internal
#' @noRd
second_der <- function(x, y) {
  deriv(middle_pts(x), deriv(x, y))
}

#' Fast Fourier transform low-pass filter
#'
#' Applies a conjugate-symmetric 4th-order Butterworth low-pass filter in the
#' frequency domain using reflection padding to prevent boundary wrap-around
#' artifacts and Gibbs ringing around sharp chromatographic peaks.
#'
#' @param x Numeric vector of signal intensities.
#' @param components Numeric(1) fraction of low-frequency components to retain
#'   (default `0.01`).
#'
#' @return Numeric vector of FFT-filtered intensities.
#' @keywords internal
#' @noRd
filter_fft <- function(x, components = 0.01) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  n <- length(x)
  if (n <= 2L || components >= 0.5) {
    return(x)
  }
  x_pad <- c(x, rev(x))
  n_pad <- length(x_pad)
  temp <- stats::fft(x_pad)
  cutoff <- max(2, n_pad * components)
  freq_idx <- pmin(seq_len(n_pad) - 1L, n_pad - seq_len(n_pad) + 1L)
  butterworth <- 1 / (1 + (freq_idx / cutoff)^8)
  out <- Re(stats::fft(temp * butterworth, inverse = TRUE)) / n_pad
  out[seq_len(n)]
}

#' Zero-phase symmetric centered running mean
#'
#' @param x Numeric vector.
#' @param k Integer(1) window width.
#'
#' @return Numeric vector of symmetrically smoothed values of length
#'   `length(x)`.
#' @keywords internal
#' @noRd
runmean_symmetric <- function(x, k) {
  k <- as.integer(k)
  if (k <= 1L || length(x) == 0L) {
    return(as.vector(x))
  }
  if (k %% 2L == 1L) {
    return(runmean_center(x, k))
  }
  0.5 * (runmean_center(x, k) + rev(runmean_center(rev(x), k)))
}

#' Signal sharpening via second and fourth derivatives
#'
#' Enhances chromatographic peak resolution via zero-phase symmetric smoothing
#' and time-aligned centered second and fourth finite differences on the
#' interior grid `3:(n - 2)`.
#'
#' @param time Numeric vector of retention times.
#' @param intensity Numeric vector of signal intensities.
#' @param k2 Numeric(1) second-derivative scaling denominator (default `250`).
#' @param k4 Numeric(1) fourth-derivative scaling denominator (default
#'   `1250000`).
#' @param sigma Numeric(1) sharpening weight factor (default `0.05`).
#' @param Smoothing_width Integer(1) running-mean window width (default `8`).
#' @param Baseline_adjust Numeric(1) additive baseline offset (default `0`).
#'
#' @return Numeric vector of sharpened intensities on the interior grid
#'   `3:(length(time) - 2)`.
#' @keywords internal
#' @noRd
signal_sharpening <- function(
  time,
  intensity,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  Smoothing_width = 8, # nolint: object_name_linter.
  Baseline_adjust = 0 # nolint: object_name_linter.
) {
  n <- length(intensity)
  if (n < 5L) {
    return(numeric(0L))
  }
  smooth_1 <- runmean_symmetric(intensity, k = Smoothing_width) +
    Baseline_adjust
  smooth_2 <- runmean_symmetric(smooth_1, k = Smoothing_width)
  deriv_2 <- second_der(x = time, y = smooth_2)
  smooth_3 <- runmean_symmetric(deriv_2, k = Smoothing_width)
  deriv_4 <- second_der(x = time[2:(n - 1L)], y = smooth_3)
  smooth_4 <- runmean_symmetric(deriv_4, k = Smoothing_width)

  sharpened <- smooth_1[3:(n - 2L)] -
    (sigma / k2 * smooth_3[2:(n - 3L)]) +
    (sigma / k4 * smooth_4)
  sharpened[is.na(sharpened)] <- 0
  sharpened
}

#' Estimate optimal retention time shift between two chromatograms via FFT
#' cross-correlation
#'
#' @param ref Reference chromatogram `data.frame` (`rtime`, `intensity`).
#' @param target Target chromatogram `data.frame` (`rtime`, `intensity`) to
#'   align to `ref`.
#' @param max_shift Numeric(1) maximum allowed absolute shift in minutes
#'   (default `0.5`).
#' @param grid_step Numeric(1) interpolation step in minutes (default `0.002`).
#'
#' @return Numeric(1) optimal retention time shift `delta` such that
#'   `target$rtime + delta` aligns with `ref$rtime`.
#' @keywords internal
#' @noRd
estimate_rt_shift <- function(
  ref,
  target,
  max_shift = 0.5,
  grid_step = 0.002
) {
  t_min <- max(min(ref$rtime, na.rm = TRUE), min(target$rtime, na.rm = TRUE))
  t_max <- min(max(ref$rtime, na.rm = TRUE), max(target$rtime, na.rm = TRUE))
  if (!is.finite(t_min) || !is.finite(t_max) || (t_max - t_min) <= grid_step) {
    return(0)
  }
  grid <- seq(t_min, t_max, by = grid_step)
  ng <- length(grid)
  if (ng < 4L) {
    return(0)
  }
  y_ref <- stats::approx(ref$rtime, ref$intensity, xout = grid, rule = 2)$y
  y_tgt <- stats::approx(
    target$rtime,
    target$intensity,
    xout = grid,
    rule = 2
  )$y
  y_ref <- y_ref - mean(y_ref, na.rm = TRUE)
  y_tgt <- y_tgt - mean(y_tgt, na.rm = TRUE)
  if (stats::sd(y_ref) == 0 || stats::sd(y_tgt) == 0) {
    return(0)
  }
  n_pad <- 2L * ng
  f_ref <- stats::fft(c(y_ref, rep(0, ng)))
  f_tgt <- stats::fft(c(y_tgt, rep(0, ng)))
  xc <- Re(stats::fft(f_ref * Conj(f_tgt), inverse = TRUE)) / n_pad
  max_lag <- min(ng - 1L, max(1L, as.integer(round(max_shift / grid_step))))
  pos_lags <- seq.int(0L, max_lag)
  neg_lags <- seq.int(-max_lag, -1L)
  lags <- c(neg_lags, pos_lags)
  xc_vals <- c(xc[n_pad + neg_lags + 1L], xc[pos_lags + 1L])
  lags[which.max(xc_vals)] * grid_step
}

#' Improve chromatogram signal via FFT filtering, resampling, and sharpening
#'
#' @param df Data frame with `rtime` and `intensity` columns.
#' @param fourier_components Numeric(1) FFT filter cutoff fraction.
#' @param frequency Numeric(1) detector sampling frequency in Hz.
#' @param resample Numeric(1) resampling multiplier.
#' @param time_min,time_max Numeric(1) retention time window bounds.
#' @param intensity_floor Numeric(1) positive shift floor for non-positive
#'   signals.
#' @param k2,k4,sigma,smoothing_width Signal sharpening parameters passed to
#'   `signal_sharpening()`.
#'
#' @return A `data.frame` with columns `rtime` and `intensity`.
#' @keywords internal
#' @noRd
improve_signal <- function(
  df,
  fourier_components = 0.01,
  frequency = 2,
  resample = 1,
  time_min = 0,
  time_max = Inf,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8
) {
  min_intensity <- min(df$intensity, na.rm = TRUE)
  shift_amount <- if (min_intensity <= 0) {
    abs(min_intensity) + intensity_floor
  } else {
    0
  }

  intensity_fourier <- filter_fft(
    x = df$intensity + shift_amount,
    components = fourier_components
  )

  f <- stats::approxfun(
    x = df$rtime,
    y = intensity_fourier
  )

  time <- seq(
    from = time_min,
    to = min(max(df$rtime), time_max),
    by = 1 / (frequency * 60 * resample)
  )

  intensity <- f(time)

  intensity_sharpened <- signal_sharpening(
    time = time,
    intensity = intensity,
    k2 = k2,
    k4 = k4,
    sigma = sigma,
    Smoothing_width = smoothing_width
  )

  n_t <- length(time)
  data.frame(
    rtime = time[3:(n_t - 2L)],
    intensity = intensity_sharpened
  )
}

#' Apply signal improvement across a list of chromatograms
#'
#' @param xs List of chromatogram data frames (`rtime`, `intensity`).
#' @param fourier_components,frequency,resample,time_min,time_max
#'   Signal processing parameters forwarded to `improve_signal()`.
#' @param intensity_floor,k2,k4,sigma,smoothing_width
#'   Additional signal processing parameters forwarded to `improve_signal()`.
#'
#' @return List of improved chromatogram data frames.
#' @keywords internal
#' @noRd
improve_signals_progress <- function(
  xs,
  fourier_components = 0.01,
  frequency = 2,
  resample = 1,
  time_min = 0,
  time_max = Inf,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8
) {
  lapply(
    xs,
    improve_signal,
    fourier_components = fourier_components,
    frequency = frequency,
    resample = resample,
    time_min = time_min,
    time_max = time_max,
    intensity_floor = intensity_floor,
    k2 = k2,
    k4 = k4,
    sigma = sigma,
    smoothing_width = smoothing_width
  )
}

#' Baseline correction of a chromatogram
#'
#' @param df Data frame with `rtime` and `intensity` columns.
#' @param method Character(1) baseline method passed to [baseline::baseline()].
#' @param ... Additional arguments passed to [baseline::baseline()].
#'
#' @return Data frame `df` with baseline-corrected `intensity`.
#' @keywords internal
#' @noRd
baseline_chromatogram <- function(df, method = "peakDetection", ...) {
  df2 <- df
  intensity <- df$intensity
  intensity[is.na(intensity)] <- 0
  intensity_baseline <- baseline::baseline(
    spectra = t(intensity),
    method = method,
    ...
  )
  df2$intensity <- as.numeric(t(intensity_baseline@corrected))
  df2
}

#' Normalize chromatogram retention time and intensity
#'
#' @param list Data frame or list with `rtime` and `intensity` elements.
#' @param shift Numeric(1) retention time shift in minutes.
#' @param normalize_intensity Logical(1); divide `intensity` by its maximum.
#' @param normalize_time Logical(1); divide `rtime` by its maximum and append
#'   `time_2`.
#'
#' @return A `data.frame` with normalized `rtime` and `intensity` (plus `time_2`
#'   when `normalize_time = TRUE`).
#' @keywords internal
#' @noRd
normalize_chromatograms_list <- function(
  list,
  shift = 0,
  normalize_intensity = TRUE,
  normalize_time = FALSE
) {
  rtime <- list$rtime + shift
  intensity <- list$intensity
  if (normalize_intensity) {
    intensity <- intensity / max(intensity)
  }
  if (normalize_time) {
    time_2 <- max(rtime)
    rtime <- rtime / time_2
    return(data.frame(rtime = rtime, intensity = intensity, time_2 = time_2))
  }
  data.frame(rtime = rtime, intensity = intensity)
}

#' Standardize chromatogram column names
#'
#' @param df Data frame containing retention time and intensity columns.
#' @param name_rt Character(1) source retention time column name.
#' @param name_intensity Character(1) source intensity column name.
#'
#' @return A `data.frame` with columns `rtime` and `intensity`.
#' @keywords internal
#' @noRd
change_intensity_name <- function(
  df,
  name_rt = "rtime",
  name_intensity = "intensity"
) {
  data.frame(
    rtime = df[[name_rt]],
    intensity = df[[name_intensity]]
  )
}

#' Extract a specific detector chromatogram from a list
#'
#' @param list Named list of detector chromatogram data frames.
#' @param type Character(1) detector key (`"bpi"`, `"cad"`, or `"pda"`).
#' @param headers Named character vector mapping detector types to header IDs.
#'
#' @return A `data.frame` with columns `rtime` and `intensity`.
#' @keywords internal
#' @noRd
extract_chromatogram <- function(list, type, headers) {
  stopifnot(
    "type must be one of 'bpi, 'cad',or 'pda'" = type %in%
      c("bpi", "cad", "pda")
  )
  name <- as.character(headers[[type]])
  idx <- which(headers == name)[1L]
  change_intensity_name(list[[idx]])
}

#' Preprocess original, improved, and baselined chromatograms
#'
#' @param detector Character(1) detector label (default `"cad"`).
#' @param fourier_components,frequency,resample,time_min,time_max
#'   Signal improvement parameters forwarded to `improve_signals_progress()`.
#' @param intensity_floor,k2,k4,sigma,smoothing_width
#'   Additional parameters forwarded to `improve_signals_progress()`.
#' @param list List of raw chromatogram data frames.
#' @param name Character vector of sample names.
#' @param shift Numeric(1) retention time shift in minutes.
#' @param baseline_method Character(1) baseline correction method.
#' @param improve_signal Logical(1); whether to apply FFT and derivative
#'   sharpening.
#'
#' @return Named list with `chromatograms_original`,
#'   `chromatograms_original_long`, `chromatograms_improved`,
#'   `chromatograms_improved_long`, `chromatograms_baselined`, and
#'   `chromatograms_baselined_long`.
#' @keywords internal
#' @noRd
preprocess_chromatograms <- function(
  detector = "cad",
  fourier_components = 0.01,
  frequency = 2,
  list,
  name,
  resample = 1,
  shift = 0,
  time_min = 0,
  time_max = Inf,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8,
  baseline_method = "peakDetection",
  improve_signal = TRUE
) {
  message("preprocessing ", detector, " chromatograms")
  chromatograms_original <- list

  if (improve_signal) {
    message("improving chromatograms")
    chromatograms_improved <- improve_signals_progress(
      xs = chromatograms_original,
      fourier_components = fourier_components,
      frequency = frequency,
      resample = resample,
      time_min = time_min,
      time_max = time_max,
      intensity_floor = intensity_floor,
      k2 = k2,
      k4 = k4,
      sigma = sigma,
      smoothing_width = smoothing_width
    )
  } else {
    message("skipping signal improvement (using original chromatograms)")
    chromatograms_improved <- lapply(chromatograms_original, function(x) {
      sel <- x$rtime >= time_min & x$rtime <= time_max
      data.frame(rtime = x$rtime[sel], intensity = x$intensity[sel])
    })
  }

  names(chromatograms_original) <- name
  names(chromatograms_improved) <- name

  chromatograms_original_long <- data.table::rbindlist(
    lapply(chromatograms_original, data.table::as.data.table),
    idcol = "id"
  )
  chromatograms_original_long[, rtime := rtime + shift]
  chromatograms_original_long[, intensity := intensity - min(intensity)]
  chromatograms_original_long[, `:=`(rt_1 = rtime, rt_2 = rtime)]

  chromatograms_improved_long <- data.table::rbindlist(
    lapply(chromatograms_improved, data.table::as.data.table),
    idcol = "id"
  )
  chromatograms_improved_long[, rtime := rtime + shift]
  chromatograms_improved_long[, `:=`(rt_1 = rtime, rt_2 = rtime)]

  message("baselining chromatograms")
  chromatograms_baselined <- lapply(
    chromatograms_improved,
    baseline_chromatogram,
    method = baseline_method
  )

  chromatograms_baselined_long <- data.table::rbindlist(
    lapply(chromatograms_baselined, data.table::as.data.table),
    idcol = "id"
  )
  chromatograms_baselined_long[,
    intensity := pmax(0, intensity - max(0, min(intensity, na.rm = TRUE)))
  ]
  chromatograms_baselined_long[, `:=`(rt_1 = rtime, rt_2 = rtime)]

  list(
    chromatograms_original = chromatograms_original,
    chromatograms_original_long = chromatograms_original_long,
    chromatograms_improved = chromatograms_improved,
    chromatograms_improved_long = chromatograms_improved_long,
    chromatograms_baselined = chromatograms_baselined,
    chromatograms_baselined_long = chromatograms_baselined_long
  )
}
