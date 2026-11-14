#' Detect peaks on a single chromatogram xy data.frame
#'
#' @param df_xy Data frame with `rtime` and `intensity` columns.
#' @param sd_max Numeric(1) maximum peak standard deviation (default `50`).
#' @param max_iter Integer(1) maximum fitting iterations (default `1000`).
#' @param noise_threshold Numeric(1) minimum peak amplitude threshold (default
#'   `0.001`).
#' @param fit Character(1) peak fitting model (`"egh"`, `"gaussian"`,
#'   or `"raw"`).
#'
#' @return A `data.table` of detected peaks with `peak_id`, `peak_max`,
#'   `rt_apex`, `rt_min`, and `rt_max`.
#' @keywords internal
#' @noRd
peaks_progress <- function(
  df_xy,
  sd_max = 50,
  max_iter = 1000,
  noise_threshold = 0.001,
  fit = "egh"
) {
  sel <- df_xy$rtime >= 0
  mat <- matrix(
    df_xy$intensity[sel],
    ncol = 1L,
    dimnames = list(as.character(df_xy$rtime[sel]), "666")
  )

  pks <- get_peaks(
    list("666" = mat),
    lambdas = "666",
    fit = fit,
    sd.max = sd_max,
    max.iter = max_iter,
    noise_threshold = noise_threshold
  )[["666"]][["666"]]

  dt <- data.table::as.data.table(pks)
  dt[, `:=`(
    peak_id = seq_len(.N),
    peak_max = height,
    rt_apex = rt,
    rt_min = start,
    rt_max = end
  )]
  dt
}

#' Compute trapezoidal numerical integral of chromatogram slice
#'
#' @param rt Numeric vector of retention times.
#' @param int Numeric vector of intensities.
#'
#' @return Numeric(1) trapezoidal integral.
#' @keywords internal
#' @noRd
trapz_integral <- function(rt, int) {
  ok <- !is.na(rt) & !is.na(int)
  rt <- rt[ok]
  int <- pmax(0, int[ok])
  n <- length(int)
  if (n == 0L) {
    return(0)
  }
  if (n == 1L) {
    return(int[1L])
  }
  ord <- order(rt)
  rt <- rt[ord]
  int <- int[ord]
  sum(diff(rt) * (int[-1L] + int[-n]) / 2)
}

#' Join detected peaks with chromatogram traces and compute peak area integral
#'
#' @param chromatograms `data.table` or `data.frame` of long chromatogram traces
#'   with `rt_1`, `rt_2`, `intensity`, and `id`.
#' @param peaks `data.table` or `data.frame` of detected peaks with `rt_min` and
#'   `rt_max`.
#' @param min_area Numeric(1) minimum relative peak area fraction to retain.
#'
#' @return A `data.table` of integrated peaks filtered by `min_area`.
#' @keywords internal
#' @noRd
join_peaks <- function(chromatograms, peaks, min_area) {
  peaks_dt <- data.table::as.data.table(peaks)
  chrom_dt <- data.table::as.data.table(chromatograms)
  data.table::setkey(peaks_dt, rt_min, rt_max)
  data.table::setkey(chrom_dt, rt_1, rt_2)

  ov <- data.table::foverlaps(peaks_dt, chrom_dt)
  ov[, integral := trapz_integral(rt_1, intensity), by = c("peak_id", "id")]
  res <- unique(
    ov[, list(peak_id, id, peak_max, rt_apex, rt_min, rt_max, integral)]
  )
  res[,
    keep_area := {
      total_int <- sum(integral, na.rm = TRUE)
      if (!is.finite(total_int) || total_int <= 0) {
        FALSE
      } else {
        (integral / total_int) >= min_area
      }
    },
    by = "id"
  ]
  res <- res[which(keep_area)][, keep_area := NULL]
  res
}

#' Normalize chromatogram slice within peak retention time bounds
#'
#' @param x Data frame or `data.table` containing `rt_min` and `rt_max`.
#' @param df_xy Full chromatogram `data.frame` (`rtime`, `intensity`).
#' @param intensity_threshold Numeric(1) minimum normalized intensity cutoff
#'   (default `0.1`).
#'
#' @return A `data.table` with `rtime` and normalized `intensity`.
#' @keywords internal
#' @noRd
normalize_chromato <- function(x, df_xy, intensity_threshold = 0.1) {
  sel <- df_xy$rtime >= x$rt_min[1L] & df_xy$rtime <= x$rt_max[1L]
  rt <- df_xy$rtime[sel]
  int <- df_xy$intensity[sel]
  rt_bounds <- if (length(rt) > 0L) c(min(rt), max(rt)) else numeric(0L)
  if (length(int) > 0L) {
    min_i <- min(int, na.rm = TRUE)
    max_i <- max(int, na.rm = TRUE)
    if (!is.finite(min_i) || !is.finite(max_i) || max_i <= min_i) {
      int <- rep(0, length(int))
    } else {
      int <- (int - min_i) / (max_i - min_i)
    }
  }
  keep <- !is.na(int) & int >= intensity_threshold
  out <- data.table::data.table(rtime = rt[keep], intensity = int[keep])
  data.table::setattr(out, "rt_bounds", rt_bounds)
  out
}

#' Prepare mzmine feature table for retention time overlap joining
#'
#' @param df Feature table `data.frame` or `data.table` in mzmine format.
#' @param min_intensity Numeric(1) minimum `intensity_max` threshold.
#' @param name Character(1) sample datafile name.
#'
#' @return A keyed `data.table` of detected features above `min_intensity`.
#' @keywords internal
#' @noRd
prepare_features <- function(df, min_intensity, name) {
  message("... preparing features")
  dt <- data.table::as.data.table(df)
  col_map <- c(
    feature_id = "id",
    rt = "rt",
    mz = "mz",
    area = "area",
    status = paste0("datafile:", name, ":feature_state"),
    rt_1 = paste0("datafile:", name, ":rt_range:min"),
    rt_2 = paste0("datafile:", name, ":rt_range:max"),
    mz_min = paste0("datafile:", name, ":mz_range:min"),
    mz_max = paste0("datafile:", name, ":mz_range:max"),
    intensity_min = paste0("datafile:", name, ":intensity_range:min"),
    intensity_max = paste0("datafile:", name, ":intensity_range:max")
  )

  sub_dt <- dt[, unname(col_map), with = FALSE]
  data.table::setnames(sub_dt, names(col_map))
  sub_dt <- sub_dt[status == "DETECTED"][, status := NULL]
  for (col in names(sub_dt)) {
    data.table::set(sub_dt, j = col, value = as.numeric(sub_dt[[col]]))
  }

  message("... keeping features above desired intensity only")
  df_features <- sub_dt[intensity_max >= min_intensity]

  message("setting joining keys")
  data.table::setkey(df_features, rt_1, rt_2)
  df_features
}

#' Normalize peak retention time to `[0, 1]` for chromatogram shape comparison
#'
#' @param x Data frame or `data.table` with `rtime` and `intensity`.
#'
#' @return Named list with `intensity` and `[0, 1]`-scaled `rtime`.
#' @keywords internal
#' @noRd
prepare_peaks <- function(x) {
  rt <- x$rtime
  if (length(rt) > 0L) {
    bounds <- attr(x, "rt_bounds")
    if (length(bounds) == 2L && bounds[2L] > bounds[1L]) {
      min_rt <- bounds[1L]
      max_rt <- bounds[2L]
    } else {
      min_rt <- min(rt)
      max_rt <- max(rt)
    }
    if (max_rt > min_rt) {
      rt <- (rt - min_rt) / (max_rt - min_rt)
    }
  }
  list(intensity = x$intensity, rtime = rt)
}

#' Prepare retention time window in seconds
#'
#' @param x Data frame or `data.table` with `rt_min` and `rt_max` in minutes.
#' @param shift Numeric(1) retention time shift in minutes.
#'
#' @return A 1-row numeric matrix with columns `rtmin` and `rtmax` in seconds.
#' @keywords internal
#' @noRd
prepare_rt <- function(x, shift = 0) {
  cbind(
    rtmin = (x$rt_min[1L] + shift) * 60,
    rtmax = (x$rt_max[1L] + shift) * 60
  )
}

#' Prepare m/z windows for a peak's features
#'
#' @param x Data frame or `data.table` with `mz_min` and `mz_max` columns.
#'
#' @return List of 1-row numeric matrices (`mzmin`, `mzmax`).
#' @keywords internal
#' @noRd
prepare_mz <- function(x) {
  n <- nrow(x)
  mz_min <- x$mz_min
  mz_max <- x$mz_max
  lapply(seq_len(n), function(z) {
    cbind(mzmin = mz_min[z], mzmax = mz_max[z])
  })
}

#' Preprocess detector peaks and overlap with MS features
#'
#' @param detector Character(1) detector name (default `"cad"`).
#' @param df_features Keyed `data.table` of MS features from
#'   `prepare_features()`.
#' @param df_long Long-format `data.table` of chromatogram points.
#' @param df_xy Single-chromatogram `data.frame` (`rtime`, `intensity`).
#' @param name Character(1) sample name.
#' @param shift Numeric(1) retention time shift in minutes.
#' @param min_area Numeric(1) minimum relative peak area fraction.
#' @param sd_max,max_iter,noise_threshold,fit Peak detection parameters passed
#'   to `peaks_progress()`.
#' @param intensity_threshold Numeric(1) normalized intensity threshold.
#'
#' @return Named list with `list_df_features_with_peaks_long`,
#'   `list_chromato_peaks`, `list_rtr`, `list_mzr`, and
#'   `df_features_without_peaks`.
#' @keywords internal
#' @noRd
preprocess_peaks <- function(
  detector = "cad",
  df_features,
  df_long,
  df_xy,
  name,
  shift = 0,
  min_area = 0,
  sd_max = 50,
  max_iter = 1000,
  noise_threshold = 0.001,
  fit = "egh",
  intensity_threshold = 0.1
) {
  message("preprocessing ", detector, " peaks")
  peaks <- peaks_progress(
    df_xy = df_xy,
    sd_max = sd_max,
    max_iter = max_iter,
    noise_threshold = noise_threshold,
    fit = fit
  )

  peaks_long <- data.table::as.data.table(peaks)
  peaks_long[, id := as.character(seq_len(.N))]
  data.table::setcolorder(peaks_long, c("id", setdiff(names(peaks_long), "id")))

  message("joining peaks")
  df_peaks <- join_peaks(
    chromatograms = df_long,
    peaks = peaks_long,
    min_area = min_area
  )

  data.table::setkey(df_peaks, rt_min, rt_max)

  message("joining within given rt tolerance")
  df_features_peaks <- data.table::foverlaps(df_features, df_peaks)

  df_features_with_peaks <- unique(
    df_features_peaks[!is.na(peak_id), !c("rt_1", "rt_2"), with = FALSE]
  )

  message("selecting features outside peaks")
  df_features_without_peaks <- unique(df_features_peaks[is.na(peak_id)])

  message("splitting by file")
  if (nrow(df_features_with_peaks) == 0L) {
    list_feat_with_pks_long <- list()
  } else {
    by_id <- split(
      df_features_with_peaks,
      factor(
        df_features_with_peaks$id,
        levels = unique(df_features_with_peaks$id)
      )
    )
    message("splitting by peak")
    by_peak <- lapply(by_id, function(x) {
      split(x, factor(x$peak_id, levels = unique(x$peak_id)))
    })
    list_feat_with_pks_long <- unname(unlist(
      by_peak,
      recursive = FALSE
    ))
  }

  if (nrow(df_features_with_peaks) == 0L) {
    message("splitting by peak")
  }

  message("normalizing chromato")
  list_chromato_with_peak <- lapply(
    list_feat_with_pks_long,
    normalize_chromato,
    df_xy = df_xy,
    intensity_threshold = intensity_threshold
  )

  message("preparing peaks chromato")
  list_chromato_peaks <- lapply(list_chromato_with_peak, prepare_peaks)

  message("preparing rt")
  list_rtr <- lapply(
    list_feat_with_pks_long,
    prepare_rt,
    shift = shift
  )

  message("preparing mz")
  list_mzr <- lapply(list_feat_with_pks_long, prepare_mz)

  list(
    list_df_features_with_peaks_long = list_feat_with_pks_long,
    list_chromato_peaks = list_chromato_peaks,
    list_rtr = list_rtr,
    list_mzr = list_mzr,
    df_features_without_peaks = df_features_without_peaks
  )
}

#' Transform extracted MS chromatogram slice to normalized rtime and intensity
#'
#' @param x List of raw feature XIC lists (`intensity`, `rtime`).
#' @param min_intensity Numeric(1) minimum normalized intensity threshold
#'   (default `0.1`).
#'
#' @return List of normalized `data.table`s (`intensity`, `rtime`).
#' @keywords internal
#' @noRd
transform_ms <- function(x, min_intensity = 0.1, preserve_window = FALSE) {
  custom_min <- function(v) if (length(v) > 0L) min(v) else Inf
  custom_max <- function(v) if (length(v) > 0L) max(v) else Inf

  lapply(seq_along(x), function(z) {
    el <- x[[z]]
    int <- el$intensity
    rt <- el$rtime
    keep <- !is.na(int)
    int <- int[keep]
    rt <- rt[keep]
    if (length(int) == 0L || custom_max(int) == custom_min(int)) {
      return(data.table::data.table(
        intensity = numeric(0L),
        rtime = numeric(0L)
      ))
    }
    if (preserve_window && custom_max(rt) > custom_min(rt)) {
      rt <- (rt - custom_min(rt)) / (custom_max(rt) - custom_min(rt))
    }
    int <- (int - custom_min(int)) / (custom_max(int) - custom_min(int))
    keep2 <- !is.na(int) & int >= min_intensity
    int <- int[keep2]
    rt <- rt[keep2]
    if (
      !preserve_window && length(rt) > 0L && custom_max(rt) > custom_min(rt)
    ) {
      rt <- (rt - custom_min(rt)) / (custom_max(rt) - custom_min(rt))
    }
    ord <- order(rt)
    data.table::data.table(intensity = int[ord], rtime = rt[ord])
  })
}

#' Extract MS1 XICs across all peaks using binary search on sorted spectra
#'
#' @param xs Integer vector of peak indices.
#' @param ms_data Spectral data object (`list(rt, mz, intensity)` or `MSnExp`).
#' @param rts List of retention time window matrices (`rtmin`, `rtmax`).
#' @param mzs List of lists of m/z window matrices (`mzmin`, `mzmax`).
#' @param nrows Optional integer(1) unused compatibility parameter.
#'
#' @return List of normalized feature XIC `data.table` lists per peak.
#' @keywords internal
#' @noRd
extract_ms_progress <- function(xs, ms_data, rts, mzs, nrows = NULL) {
  spectra <- parse_ms_spectra(ms_data)
  sp_rt <- spectra$rt
  sp_mz <- spectra$mz
  sp_int <- spectra$intensity

  lapply(xs, function(x) {
    message("CAD Peak: ", x)
    rt_mat <- rts[[x]]
    mz_mat <- do.call(rbind, mzs[[x]])
    n_feat <- nrow(mz_mat)
    mz1 <- mz_mat[, 1L]
    mz2 <- mz_mat[, 2L]

    # Binary search on sorted retention times
    r_rt <- findInterval(rt_mat[1L, 2L], sp_rt)
    l_rt <- findInterval(rt_mat[1L, 1L], sp_rt, left.open = TRUE)
    if (r_rt <= l_rt) {
      raw_xics <- replicate(
        n_feat,
        list(intensity = numeric(0L), rtime = numeric(0L)),
        simplify = FALSE
      )
      return(transform_ms(raw_xics, preserve_window = TRUE))
    }

    scan_idx <- (l_rt + 1L):r_rt
    rt_slice <- sp_rt[scan_idx]
    n_scans <- length(scan_idx)

    # Extract intensity matrix (n_feat x n_scans) via binary search per scan
    int_mat <- matrix(NA_real_, nrow = n_feat, ncol = n_scans)
    for (s in seq_len(n_scans)) {
      j <- scan_idx[s]
      m <- sp_mz[[j]]
      if (length(m) == 0L) {
        next
      }
      int_mat[, s] <- 0
      v <- sp_int[[j]]
      r_mz <- findInterval(mz2, m)
      l_mz <- findInterval(mz1, m, left.open = TRUE)
      n_hit <- r_mz - l_mz
      one <- which(n_hit == 1L)
      if (length(one) > 0L) {
        int_mat[one, s] <- v[r_mz[one]]
      }
      many <- which(n_hit > 1L)
      if (length(many) > 0L) {
        for (k in many) {
          int_mat[k, s] <- sum(v[(l_mz[k] + 1L):r_mz[k]])
        }
      }
    }

    raw_xics <- lapply(seq_len(n_feat), function(f_idx) {
      list(intensity = int_mat[f_idx, ], rtime = rt_slice)
    })
    transform_ms(raw_xics, preserve_window = TRUE)
  })
}

#' Prepare extracted MS peaks for correlation comparison
#'
#' @param x List of normalized feature XIC `data.table`s.
#'
#' @return List of `list(intensity, rtime)` pairs.
#' @keywords internal
#' @noRd
extract_ms_peak <- function(x) {
  lapply(seq_along(x), function(z) {
    list(intensity = x[[z]]$intensity, rtime = x[[z]]$rtime)
  })
}

#' Compare CAD peak shape with MS feature XIC peak shapes
#'
#' @param x Integer(1) peak index.
#' @param list_ms_peaks List of extracted MS peak lists.
#' @param peaks_prelist Output list from `preprocess_peaks()`.
#'
#' @return List of Pearson correlation coefficients for each feature in
#'   peak `x`.
#' @keywords internal
#' @noRd
compare_peaks <- function(x, list_ms_peaks, peaks_prelist) {
  ms_pks <- list_ms_peaks[[x]]
  if (length(ms_pks) == 0L) {
    return(rep(list(0), length(x)))
  }
  chrom_pk <- peaks_prelist$list_chromato_peaks[[x]]
  x_rt <- chrom_pk$rtime
  x_int <- chrom_pk$intensity

  safe_cor <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 2L) {
      return(0)
    }
    a_ok <- a[ok]
    b_ok <- b[ok]
    if (stats::sd(a_ok) == 0 || stats::sd(b_ok) == 0) {
      return(0)
    }
    r_val <- stats::cor(a_ok, b_ok)
    if (is.na(r_val)) 0 else r_val
  }

  lapply(seq_along(ms_pks), function(z) {
    y_pk <- ms_pks[[z]]
    y_rt <- y_pk$rtime
    y_int <- y_pk$intensity
    if (length(y_int) <= 1L || length(x_int) <= 1L) {
      return(0)
    }
    if (length(x_rt) != length(y_rt) || !all(x_rt == y_rt)) {
      tol <- min(mean(diff(x_rt)), mean(diff(y_rt)))
      x_min <- min(x_rt)
      x_max <- max(x_rt)
      y_min <- min(y_rt)
      y_max <- max(y_rt)
      if (abs(x_min - y_min) <= tol && abs(x_max - y_max) <= tol) {
        aligned_int <- stats::approx(x_rt, x_int, xout = y_rt, rule = 2)$y
        return(safe_cor(aligned_int, y_int))
      }
      grid_rt <- sort(unique(c(x_rt, y_rt)))
      x_eval <- stats::approx(
        x_rt,
        x_int,
        xout = grid_rt,
        yleft = 0,
        yright = 0
      )$y
      y_eval <- stats::approx(
        y_rt,
        y_int,
        xout = grid_rt,
        yleft = 0,
        yright = 0
      )$y
      return(safe_cor(x_eval, y_eval))
    }
    safe_cor(x_int, y_int)
  })
}
