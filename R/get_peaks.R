#' Detect and integrate chromatographic peaks
#'
#' @param chrom_list Matrix or list of chromatogram intensity matrices (rows =
#'   timepoints, columns = wavelengths/channels).
#' @param lambdas Character or numeric vector of wavelength/channel identifiers.
#' @param fit Character(1) peak fitting model: `"egh"` (exponential-Gaussian
#'   hybrid), `"gaussian"`, or `"raw"`.
#' @param sd.max Numeric(1) maximum allowed peak width (`sd`) in index units.
#' @param max.iter Integer(1) maximum fitting iterations.
#' @param time.units Character(1) output time unit scaling (`"min"`, `"s"`, or
#'   `"ms"`).
#' @param estimate_purity Logical(1) compatibility parameter.
#' @param noise_threshold Numeric(1) minimum peak amplitude threshold.
#' @param collapse Logical(1); bind wavelength tables per sample into a single
#'   `data.frame`.
#' @param ... Additional arguments passed to internal `find_peaks()`.
#'
#' @return A `peak_list` object (nested list of peak `data.frame`s per sample
#'   and wavelength).
#' @author Ethan Bass
#' @note Adapted from `{chromatographR}` with parallelization and external dependencies removed.
#' @source https://github.com/ethanbass/chromatographR
#' @keywords internal
#' @noRd
get_peaks <- function(
  chrom_list,
  lambdas,
  fit = c("egh", "gaussian", "raw"),
  sd.max = 50,
  max.iter = 100,
  time.units = c("min", "s", "ms"),
  estimate_purity = FALSE,
  noise_threshold = 0.001,
  collapse = FALSE,
  ...
) {
  remove_bad_peaks <- function(pks, n) {
    pks[
      which(
        stats::complete.cases(pks[, c("rt", "start", "end")]) &
          pks[, "rt"] >= 1 &
          pks[, "rt"] <= n
      ),
      ,
      drop = FALSE
    ]
  }

  convert_indices_to_times <- function(x, chrom_list, idx, tfac) {
    get_times <- function(obj, idx = 1) {
      if (inherits(obj, "chrom_list") || inherits(obj, "list")) {
        as.numeric(rownames(obj[[idx]]))
      } else if (inherits(obj, "matrix")) {
        as.numeric(rownames(obj))
      }
    }
    timepoints <- get_times(chrom_list, idx = idx)
    tdiff <- signif(stats::median(diff(timepoints)))
    x[, c("rt", "start", "end")] <- cbind(
      timepoints[x[, "rt"]],
      timepoints[x[, "start"]],
      timepoints[x[, "end"]]
    )
    x[, c("sd", "FWHM", "area")] <- x[, c("sd", "FWHM", "area")] * tdiff * tfac
    if (!is.null(x$tau)) {
      x[, c("tau")] <- x[, c("tau")] * tdiff * tfac
    }
    x
  }

  time.units <- match.arg(time.units, c("min", "s", "ms"))
  tfac <- switch(time.units, "min" = 1, "s" = 60, "ms" = 60 * 1000)
  fit <- match.arg(fit, c("egh", "gaussian", "raw"))
  chrom_list_string <- deparse(substitute(chrom_list))
  if (is.matrix(chrom_list)) {
    chrom_list <- list(chrom_list)
  }
  if (missing(lambdas)) {
    if (ncol(chrom_list[[1L]]) == 1L) {
      lambdas <- colnames(chrom_list[[1L]])
    } else {
      stop("Wavelengths (`lambdas`) must be provided.", call. = FALSE)
    }
  }
  if (is.numeric(lambdas)) {
    lambdas <- as.character(lambdas)
  }
  if (is.null(names(chrom_list))) {
    warning(
      "Sample names not found. It is recommended to include names for your samples.",
      immediate. = TRUE
    )
    names(chrom_list) <- seq_along(chrom_list)
  }

  find_peaks <- function(
    y,
    smooth_type = c("gaussian", "box", "savgol", "mva", "tmva", "none"),
    smooth_window = 0.001,
    slope_thresh = 0,
    amp_thresh = 0,
    bounds = TRUE
  ) {
    if (!is.vector(y)) {
      stop("Please provide a vector to argument `y` to proceed.", call. = FALSE)
    }
    smooth_type <- match.arg(
      smooth_type,
      c("gaussian", "box", "savgol", "mva", "tmva", "none")
    )
    if (smooth_window < 1) {
      smooth_window <- max(length(y) * smooth_window, 1)
    }
    if (smooth_type == "savgol") {
      if ((smooth_window %% 2) == 0) {
        smooth_window <- smooth_window + 1
      }
      fc <- (smooth_window - 1) / 2
      X <- outer(-fc:fc, 0:4, FUN = "^")
      sv <- svd(X)
      pos_d <- sv$d > max(dim(X)) * .Machine$double.eps * sv$d[1L]
      Y <- sv$v[, pos_d, drop = FALSE] %*%
        ((1 / sv$d[pos_d]) * t(sv$u[, pos_d, drop = FALSE]))
      smoothed <- stats::convolve(diff(y), rev(Y[1L, ]), type = "o")
      d <- smoothed[(fc + 1):(length(smoothed) - fc)]
    } else if (smooth_type == "mva") {
      d <- runmean_center(diff(y), k = smooth_window)
    } else if (smooth_type == "gaussian") {
      d <- diff(
        stats::ksmooth(
          seq_along(y),
          y,
          kernel = "normal",
          bandwidth = smooth_window
        )$y
      )
    } else if (smooth_type == "box") {
      d <- diff(
        stats::ksmooth(
          seq_along(y),
          y,
          kernel = "box",
          bandwidth = smooth_window
        )$y
      )
    } else if (smooth_type == "tmva") {
      d <- runmean_center(
        runmean_center(diff(y), k = smooth_window),
        k = smooth_window
      )
    } else {
      d <- diff(y)
    }

    nd <- length(d)
    p1 <- which(sign(d[1:(nd - 1L)]) > sign(d[2:nd]))
    p2 <- which(abs(diff(d)) > slope_thresh)
    p3 <- which(y > amp_thresh)
    p <- intersect(intersect(p1, p2), p3)
    if (bounds) {
      p4 <- which(sign(d[1:(nd - 1L)]) < sign(d[2:nd]))
      if (length(p) == 0L || length(p4) == 0L) {
        return(data.frame(
          pos = p,
          lower = rep(1, length(p)),
          upper = rep(length(y), length(p))
        ))
      }
      idx_l <- findInterval(p - 1L, p4)
      bl <- ifelse(idx_l > 0L, p4[pmax(1L, idx_l)], 1)
      idx_u <- findInterval(p, p4) + 1L
      bu <- ifelse(
        idx_u <= length(p4),
        p4[pmin(length(p4), idx_u)],
        length(y)
      )
      p <- data.frame(pos = p, lower = bl, upper = bu)
    }
    p
  }

  gaussian_fn <- function(x, center = 0, width = 1, height = NULL, floor = 0) {
    twoVar <- 2 * width * width
    sqrt2piVar <- sqrt(pi * twoVar)
    y <- exp(-(x - center)^2 / twoVar) / sqrt2piVar
    if (!is.null(height)) {
      y <- y * sqrt2piVar * height
    }
    if (is.null(floor)) y else y + floor
  }

  egh_fn <- function(x, center, width, height, tau = 0) {
    result <- rep(0, length(x))
    index <- which(2 * width^2 + tau * (x - center) > 0)
    result[index] <- height *
      exp(-(x[index] - center)^2 / (2 * width^2 + tau * (x[index] - center)))
    result
  }

  result <- lapply(seq_along(chrom_list), function(sample) {
    suppressWarnings({
      ptable <- lapply(lambdas, function(lambda) {
        x <- chrom_list[[sample]]
        col_nums <- as.numeric(colnames(x))
        lambda.idx <- if (lambda == "max") {
          which.max(x)
        } else if (length(col_nums) == 1L) {
          1L
        } else {
          which(col_nums == as.numeric(lambda))
        }
        if (length(lambda.idx) == 0L || is.na(lambda.idx)) {
          stop(
            "The specified wavelength (`lambda`) could not be found!",
            call. = FALSE
          )
        }

        y <- x[, lambda.idx]
        pos <- find_peaks(y, ...)

        tabnames <- switch(
          fit,
          "gaussian" = c(
            "rt",
            "start",
            "end",
            "sd",
            "FWHM",
            "height",
            "area",
            "r.squared",
            "purity"
          ),
          "egh" = c(
            "rt",
            "start",
            "end",
            "sd",
            "tau",
            "FWHM",
            "height",
            "area",
            "r.squared",
            "purity"
          ),
          "raw" = c(
            "rt",
            "start",
            "end",
            "sd",
            "FWHM",
            "height",
            "area",
            "purity"
          )
        )
        noPeaksMat <- matrix(
          rep(NA_real_, length(tabnames)),
          nrow = 1L,
          dimnames = list(NULL, tabnames)
        )
        ny <- length(y)
        on.edge <- pos$pos <= 1L |
          pos$pos >= ny |
          is.na(y[pmin(ny, pos$pos + 1L)]) |
          is.na(y[pmax(1L, pos$pos - 1L)])
        pos <- pos[!on.edge, , drop = FALSE]

        if (nrow(pos) == 0L) {
          pks <- noPeaksMat
        } else {
          fitpk <- function(row_pos) {
            xloc <- row_pos[1L]
            peak.loc <- seq.int(row_pos[2L], row_pos[3L])
            y_pk <- y[peak.loc]
            if (fit == "raw") {
              area <- sum(
                diff(peak.loc) * mean(c(y_pk[-1L], utils::tail(y_pk, -1L)))
              )
              return(c(
                rt = xloc,
                start = row_pos[2L],
                end = row_pos[3L],
                sd = row_pos[3L] - row_pos[2L],
                FWHM = 2.35 * row_pos[3L] - row_pos[2L],
                height = y[xloc],
                area = area,
                purity = NA_real_
              ))
            }
            start.center <- xloc
            start.height <- y[xloc]
            start.width <- sum(y_pk > (start.height / 2)) / 2
            y_fit <- if (fit == "gaussian") {
              gaussian_fn(
                peak.loc,
                start.center,
                start.width,
                start.height,
                NULL
              )
            } else {
              egh_fn(peak.loc, start.center, start.width, start.height, 0)
            }
            area <- sum(
              diff(peak.loc) * mean(c(y_fit[-1L], utils::tail(y_fit, -1L)))
            )
            r_sq <- tryCatch(
              summary(stats::lm(y_fit ~ y_pk))$r.squared,
              error = function(e) NA_real_
            )
            if (fit == "gaussian") {
              c(
                rt = start.center,
                start = row_pos[2L],
                end = row_pos[3L],
                sd = start.width,
                FWHM = 2.35 * start.width,
                height = start.height,
                area = area,
                r.squared = r_sq,
                purity = NA_real_
              )
            } else {
              c(
                rt = start.center,
                start = row_pos[2L],
                end = row_pos[3L],
                sd = start.width,
                tau = 0,
                FWHM = 2.35 * start.width,
                height = start.height,
                area = area,
                r.squared = r_sq,
                purity = NA_real_
              )
            }
          }

          huhn <- data.frame(t(apply(pos, 1L, fitpk)))
          colnames(huhn) <- tabnames
          huhn <- data.frame(lapply(huhn, as.numeric))
          if (!is.null(sd.max)) {
            huhn <- huhn[huhn$sd < sd.max, , drop = FALSE]
          }
          pks <- try(huhn[huhn$rt > 0, , drop = FALSE], silent = TRUE)
          if (inherits(pks, "try-error")) {
            pks <- noPeaksMat
          }
        }
        pks <- cbind(sample = names(chrom_list)[sample], lambda, pks)
        pks <- remove_bad_peaks(pks, n = nrow(chrom_list[[sample]]))
        convert_indices_to_times(
          pks,
          chrom_list = chrom_list,
          idx = sample,
          tfac = tfac
        )
      })
    })
    names(ptable) <- lambdas
    if (collapse) {
      ptable <- do.call(rbind, ptable)
      rownames(ptable) <- NULL
    }
    ptable
  })
  names(result) <- names(chrom_list)
  structure(
    result,
    chrom_list = chrom_list_string,
    lambdas = lambdas,
    fit = fit,
    sd.max = sd.max,
    max.iter = max.iter,
    time.units = time.units,
    class = c("peak_list", "list")
  )
}
