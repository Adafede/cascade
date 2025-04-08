#' Get peaks
#'
#' @author Ethan Bass
#'
#' @note
#' This was imported from \{chromatographR\} package and parallelization
#' was removed as it was causing issues on Windows.
#'
#' @source https://github.com/ethanbass/chromatographR
#'
#' @param chrom_list Chrom list
#' @param lambdas Lambdas
#' @param fit Fit
#' @param sd.max Sd max
#' @param max.iter Max iter
#' @param time.units Time units
#' @param estimate_purity Estimate purity
#' @param noise_threshold Noise Threshold
#' @param collapse Collapse
#' @param ... ...
#'
#' @return Peaks
#'
#' @examples NULL
get_peaks <- function(
  chrom_list,
  lambdas,
  fit = c("egh", "gaussian", "raw"),
  sd.max = 50,
  max.iter = 100,
  time.units = c("min", "s", "ms"),
  estimate_purity = FALSE,
  noise_threshold = .001,
  collapse = FALSE,
  ...
) {
  remove_bad_peaks <- function(pks) {
    pks[
      which(
        apply(pks, 1, function(x) {
          !all(is.na(x))
        }) &
          length(pks[, "rt"] > pks[, "start"]) &
          length(pks[, "rt"] < pks[, "end"])
      ),
      ,
      drop = FALSE
    ]
  }

  convert_indices_to_times <- function(x, chrom_list, idx, tfac) {
    get_times <- function(x, idx = 1) {
      if (inherits(x, "peak_table")) {
        x <- get_chrom_list(x)
      }
      if (inherits(x, "chrom_list") | inherits(x, "list")) {
        as.numeric(rownames(x[[idx]]))
      } else if (inherits(x, "matrix")) {
        as.numeric(rownames(x))
      }
    }
    get_time_resolution <- function(chrom_list, idx = 1) {
      ts <- get_times(x = chrom_list, idx = idx)
      signif(stats::median(diff(ts)))
    }

    timepoints <- get_times(chrom_list, idx = idx)
    tdiff <- get_time_resolution(chrom_list, idx = idx)
    x[, c("rt", "start", "end")] <- sapply(
      c("rt", "start", "end"),
      function(j) {
        timepoints[x[, j]]
      }
    )
    x[, c("sd", "FWHM", "area")] <- x[, c("sd", "FWHM", "area")] * tdiff * tfac
    if (!is.null(x$tau)) {
      x[, c("tau")] <- x[, c("tau")] * tdiff * tfac
    }
    x
  }

  get_purity <- function(
    x,
    pos,
    weight = 1,
    cutoff = 0.05,
    noise_variance = NULL,
    noise_threshold = 0.01,
    lambdas,
    try = TRUE
  ) {
    get_purity_values <- function(
      x,
      pos,
      weight = 1,
      noise_variance = NULL,
      noise_threshold = 0.005,
      lambdas
    ) {
      if (missing(lambdas)) {
        lambdas <- seq_len(ncol(x))
      }
      ((1 - get_spectral_similarity(x, pos))) /
        (1 -
          get_agilent_threshold(
            x,
            pos,
            weight = weight,
            noise_variance = noise_variance,
            noise_threshold = noise_threshold,
            lambdas = lambdas
          ))
    }
    if (try) {
      try(
        {
          if (missing(lambdas)) {
            lambdas <- seq_len(ncol(x))
          }
          if (is.character(lambdas)) {
            lambdas <- which(as.numeric(colnames(x)) %in% lambdas)
          }
          p <- get_purity_values(
            x,
            pos,
            weight = weight,
            noise_variance = noise_variance,
            lambdas = lambdas
          )
          mean(p[trim_peak(x, pos, cutoff = cutoff)] < 1, na.rm = TRUE)
        },
        NA
      )
    } else {
      NA
    }
  }

  time.units <- match.arg(time.units, c("min", "s", "ms"))
  tfac <- switch(time.units, "min" = 1, "s" = 60, "ms" = 60 * 1000)
  fit <- match.arg(fit, c("egh", "gaussian", "raw"))
  chrom_list_string <- deparse(substitute(chrom_list))
  if (class(chrom_list)[1] == "matrix") {
    chrom_list <- list(chrom_list)
  }
  if (missing(lambdas)) {
    if (ncol(chrom_list[[1]]) == 1) {
      lambdas <- colnames(chrom_list[[1]])
    } else {
      stop("Wavelengths (`lambdas`) must be provided.")
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
  peaks <- list()

  result <- purrr::map(seq_along(chrom_list), function(sample) {
    suppressWarnings(
      ptable <- purrr::map(lambdas, function(lambda) {
        fit_peaks <- function(
          x,
          lambda,
          pos = NULL,
          sd.max = 50,
          fit = c("egh", "gaussian", "raw"),
          max.iter = 1000,
          estimate_purity = TRUE,
          noise_threshold = .001,
          ...
        ) {
          get_lambda_idx <- function(lambda, lambdas, y, allow_max = TRUE) {
            if (lambda == "max") {
              if (allow_max) {
                lambda.idx <- which.max(y)
              } else {
                stop(
                  "Wavelength (`lambda`) must be specified for interactive scanning."
                )
              }
            } else {
              lambda.idx <- ifelse(
                (length(lambdas) == 1),
                1,
                which(lambdas == as.numeric(lambda))
              )
            }
            if (is.na(lambda.idx) | length(lambda.idx) == 0) {
              stop("The specified wavelength (`lambda`) could not be found!")
            }
            lambda.idx
          }

          lambda.idx <- get_lambda_idx(lambda, as.numeric(colnames(x)))
          y <- x[, lambda.idx]
          fit <- match.arg(fit, c("egh", "gaussian", "raw"))
          if (is.null(pos)) {
            find_peaks <- function(
              y,
              smooth_type = c(
                "gaussian",
                "box",
                "savgol",
                "mva",
                "tmva",
                "none"
              ),
              smooth_window = .001,
              slope_thresh = 0,
              amp_thresh = 0,
              bounds = TRUE
            ) {
              if (!is.vector(y)) {
                stop("Please provide a vector to argument `y` to proceed.")
              }
              smooth_type <- match.arg(
                smooth_type,
                c("gaussian", "box", "savgol", "mva", "tmva", "none")
              )
              if (smooth_window < 1) {
                smooth_window <- max(length(y) * smooth_window, 1)
              }
              # compute derivative (with or without smoothing)
              if (smooth_type == "savgol") {
                if ((smooth_window %% 2) == 0) {
                  smooth_window <- smooth_window + 1
                }
                savgol <- function(T, fl, forder = 4, dorder = 0) {
                  stopifnot(is.numeric(T), is.numeric(fl))
                  if (fl <= 1 || fl %% 2 == 0) {
                    stop("Argument 'fl' must be an odd integer greater than 1.")
                  }
                  n <- length(T)

                  # -- calculate filter coefficients --
                  fc <- (fl - 1) / 2 # index: window left and right
                  X <- outer(-fc:fc, 0:forder, FUN = "^") # polynomial terms and coeffs
                  Y <- pinv(X) # pseudoinverse

                  # -- filter via convolution and take care of the end points --
                  T2 <- stats::convolve(T, rev(Y[(dorder + 1), ]), type = "o") # convolve(...)
                  T2 <- T2[(fc + 1):(length(T2) - fc)]

                  Tsg <- (-1)^dorder * T2
                  return(Tsg)
                }
                d <- savgol(diff(y), fl = smooth_window)
              } else if (smooth_type == "mva") {
                d <- caTools::runmean(diff(y), k = smooth_window)
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
                  ksmooth(
                    seq_along(y),
                    y,
                    kernel = "box",
                    bandwidth = smooth_window
                  )$y
                )
              } else if (smooth_type == "tmva") {
                d <- caTools::runmean(
                  caTools::runmean(diff(y), k = smooth_window),
                  k = smooth_window
                )
              } else {
                d <- diff(y)
              }

              # detect zero-crossing of first derivative (peak apex)
              p1 <- which(sign(d[1:(length(d) - 1)]) > sign(d[2:length(d)]))

              # detect second derivative exceeding slope threshold
              p2 <- which(abs(diff(d)) > slope_thresh)

              # detect y-vals exceeding amplitude threshold
              p3 <- which(y > amp_thresh)
              p <- intersect(intersect(p1, p2), p3)
              if (bounds) {
                p4 <- which(sign(d[1:(length(d) - 1)]) < sign(d[2:length(d)]))

                # find lower bound
                suppressWarnings(
                  bl <- sapply(p, function(v) {
                    max(p4[p4 < v])
                  })
                )
                bl[which(bl == -Inf)] <- 1

                # find upper bound
                suppressWarnings(
                  bu <- sapply(p, function(v) {
                    min(p4[p4 > v])
                  })
                )
                bu[which(bu == Inf)] <- length(y)
                p <- data.frame(
                  pos = p,
                  lower = bl,
                  upper = bu
                )
              }
              p
            }
            pos <- find_peaks(y, ...)
          }
          if (ncol(x) == 1) {
            estimate_purity <- FALSE
          }
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
              "r-squared",
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
            rep(NA, length(tabnames)),
            nrow = 1,
            dimnames = list(NULL, tabnames)
          )
          on.edge <- sapply(pos$pos, function(x) {
            x <= 1 || is.na(y[x + 1]) || is.na(y[x - 1])
          })
          pos <- pos[!on.edge, ]

          if (nrow(pos) == 0) {
            return(noPeaksMat)
          }

          fitpk_gaussian <- function(
            x,
            pos,
            lambda,
            max.iter,
            estimate_purity = TRUE,
            noise_threshold = .001,
            ...
          ) {
            fit_gaussian <- function(
              x,
              y,
              start.center = NULL,
              start.width = NULL,
              start.height = NULL,
              start.floor = NULL,
              fit.floor = FALSE,
              max.iter = 1000
            ) {
              gaussian <- function(
                x,
                center = 0,
                width = 1,
                height = NULL,
                floor = 0
              ) {
                # adapted from Earl F. Glynn;  Stowers Institute for Medical Research, 2007
                twoVar <- 2 * width * width
                sqrt2piVar <- sqrt(pi * twoVar)
                y <- exp(-(x - center)^2 / twoVar) / sqrt2piVar

                # by default, the height is such that the curve has unit volume
                if (!is.null(height)) {
                  scalefactor <- sqrt2piVar
                  y <- y * scalefactor * height
                }
                y + floor
              }
              # estimate starting values
              who.max <- which.max(y)
              if (is.null(start.center)) {
                start.center <- x[who.max]
              }
              if (is.null(start.height)) {
                start.height <- y[who.max]
              }
              if (is.null(start.width)) {
                start.width <- sum(y > (start.height / 2)) / 2
              }

              # call the Nonlinear Least Squares, either fitting the floor too or not
              controlList <- stats::nls.control(
                maxiter = max.iter,
                minFactor = 1 / 512,
                warnOnly = TRUE
              )
              starts <- list(
                "center" = start.center,
                "width" = start.width,
                "height" = start.height
              )
              if (!fit.floor) {
                nlsAns <- try(
                  nlsLM(
                    y ~ gaussian(x, center, width, height),
                    start = starts,
                    control = controlList
                  ),
                  silent = TRUE
                )
              } else {
                if (is.null(start.floor)) {
                  start.floor <- stats::quantile(y, seq(0, 1, 0.1))[2]
                }
                starts <- c(starts, "floor" = start.floor)
                nlsAns <- try(
                  nlsLM(
                    y ~ gaussian(x, center, width, height, floor),
                    start = starts,
                    control = controlList
                  ),
                  silent = TRUE
                )
              }

              # package up the results to pass back

              if (inherits(nlsAns, "try-error")) {
                yAns <- gaussian(
                  x,
                  start.center,
                  start.width,
                  start.height,
                  start.floor
                )
                out <- list(
                  "center" = start.center,
                  "width" = start.width,
                  "height" = start.height,
                  "y" = yAns,
                  "residual" = y - yAns
                )
                floorAns <- if (fit.floor) {
                  start.floor
                } else {
                  0
                }
              } else {
                coefs <- stats::coef(nlsAns)
                out <- list(
                  "center" = coefs[1],
                  "width" = coefs[2],
                  "height" = coefs[3],
                  "y" = stats::fitted(nlsAns),
                  "residual" = stats::residuals(nlsAns)
                )
                floorAns <- if (fit.floor) {
                  coefs[4]
                } else {
                  0
                }
              }
              if (fit.floor) {
                out <- c(out, "floor" = floorAns)
              }
              return(out)
            }
            y <- x[, lambda]
            xloc <- pos[1]
            peak.loc <- seq.int(pos[2], pos[3])
            suppressWarnings(
              m <- fit_gaussian(
                peak.loc,
                y[peak.loc],
                start.center = xloc,
                start.height = y[xloc],
                max.iter = max.iter
              )
            )
            area <- sum(diff(peak.loc) * mean(c(m$y[-1], utils::tail(m$y, -1)))) # trapezoidal integration
            r.squared <- try(
              summary(stats::lm(m$y ~ y[peak.loc]))$r.squared,
              silent = TRUE
            )
            purity <- get_purity(
              x = x,
              pos = pos,
              try = estimate_purity,
              noise_threshold = noise_threshold
            )
            c(
              "rt" = m$center,
              "start" = pos[2],
              "end" = pos[3],
              "sd" = m$width,
              "FWHM" = 2.35 * m$width,
              "height" = y[xloc],
              "area" = area,
              "r.squared" = r.squared,
              purity = purity
            )
          }

          fitpk_egh <- function(
            x,
            pos,
            lambda,
            max.iter,
            estimate_purity = TRUE,
            noise_threshold = .001
          ) {
            fit_egh <- function(
              x1,
              y1,
              start.center = NULL,
              start.width = NULL,
              start.tau = NULL,
              start.height = NULL,
              start.floor = NULL,
              fit.floor = FALSE,
              max.iter = 1000
            ) {
              egh <- function(x, center, width, height, tau, floor = 0) {
                result <- rep(0, length(x))
                index <- which(2 * width^2 + tau * (x - center) > 0)
                result[index] <- height *
                  exp(
                    -(x[index] - center)^2 /
                      (2 * width^2 + tau * (x[index] - center))
                  )
                return(result)
              }
              # try to find the best egh to fit the given data

              # make some rough estimates from the values of Y
              who.max <- which.max(y1)
              if (is.null(start.center)) {
                start.center <- x1[who.max]
              }
              if (is.null(start.height)) {
                start.height <- y1[who.max]
              }
              if (is.null(start.width)) {
                start.width <- sum(y1 > (start.height / 2)) / 2
              }
              if (is.null(start.tau)) {
                start.tau <- 0
              }
              # call the Nonlinear Least Squares, either fitting the floor too or not
              controlList <- stats::nls.control(
                maxiter = max.iter,
                minFactor = 1 / 512,
                warnOnly = TRUE
              )
              starts <- list(
                "center" = start.center,
                "width" = start.width,
                "height" = start.height,
                "tau" = start.tau
              )
              if (!fit.floor) {
                nlsAns <- try(
                  nlsLM(
                    y1 ~ egh(x1, center, width, height, tau),
                    start = starts,
                    control = controlList
                  ),
                  silent = TRUE
                )
              } else {
                if (is.null(start.floor)) {
                  start.floor <- stats::quantile(y1, seq(0, 1, 0.1))[2]
                }
                starts <- c(starts, "floor" = start.floor)
                nlsAns <- try(
                  nlsLM(
                    y1 ~ egh(x1, center, width, height, tau, floor),
                    start = starts,
                    control = controlList
                  ),
                  silent = TRUE
                )
              }

              # package up the results to pass back
              if (inherits(nlsAns, "try-error")) {
                yAns <- egh(
                  x1,
                  start.center,
                  start.width,
                  start.height,
                  start.tau,
                  start.floor
                )
                out <- list(
                  "center" = start.center,
                  "width" = start.width,
                  "height" = start.height,
                  "tau" = start.tau,
                  "y" = yAns,
                  "residual" = y1 - yAns
                )
                floorAns <- if (fit.floor) {
                  start.floor
                } else {
                  0
                }
              } else {
                coefs <- stats::coef(nlsAns)
                out <- list(
                  "center" = coefs[1],
                  "width" = coefs[2],
                  "height" = coefs[3],
                  "tau" = coefs[4],
                  "y" = stats::fitted(nlsAns),
                  "residual" = stats::residuals(nlsAns)
                )
                floorAns <- if (fit.floor) {
                  coefs[5]
                } else {
                  0
                }
              }

              if (fit.floor) {
                out <- c(out, "floor" = floorAns)
              }
              return(out)
            }
            y <- x[, lambda]
            xloc <- pos[1]
            peak.loc <- seq.int(pos[2], pos[3])
            suppressWarnings(
              m <- fit_egh(
                peak.loc,
                y[peak.loc],
                start.center = xloc,
                start.height = y[xloc],
                max.iter = max.iter
              )
            )
            r.squared <- try(
              summary(stats::lm(m$y ~ y[peak.loc]))$r.squared,
              silent = TRUE
            )
            purity <- get_purity(
              x = x,
              pos = pos,
              try = estimate_purity,
              noise_threshold = noise_threshold
            )
            # trapezoidal integration
            area <- sum(diff(peak.loc) * mean(c(m$y[-1], utils::tail(m$y, -1))))
            c(
              "rt" = m$center,
              "start" = pos[2],
              "end" = pos[3],
              "sd" = m$width,
              "tau" = m$tau,
              "FWHM" = 2.35 * m$width,
              "height" = y[xloc],
              "area" = area,
              "r.squared" = r.squared,
              purity = purity
            )
          }

          fitpk_raw <- function(
            x,
            pos,
            lambda,
            max.iter,
            estimate_purity = TRUE,
            noise_threshold = .001
          ) {
            y <- x[, lambda]
            xloc <- pos[1]
            peak.loc <- seq.int(pos[2], pos[3])

            # perform trapezoidal integration on raw signal
            area <- sum(
              diff(peak.loc) *
                mean(c(y[peak.loc][-1], utils::tail(y[peak.loc], -1)))
            )
            purity <- get_purity(
              x = x,
              pos = pos,
              try = estimate_purity,
              noise_threshold = noise_threshold
            )
            c(
              "rt" = pos[1],
              "start" = pos[2],
              "end" = pos[3],
              "sd" = pos[3] - pos[2],
              "FWHM" = 2.35 * pos[3] - pos[2],
              "height" = y[xloc],
              "area" = area,
              purity = purity
            )
          }

          fitpk <- switch(
            fit,
            "gaussian" = fitpk_gaussian,
            "egh" = fitpk_egh,
            "raw" = fitpk_raw
          )

          huhn <- data.frame(t(
            apply(
              pos,
              1,
              fitpk,
              x = x,
              lambda = lambda.idx,
              max.iter = max.iter,
              estimate_purity = estimate_purity,
              noise_threshold = noise_threshold
            )
          ))
          colnames(huhn) <- tabnames
          huhn <- data.frame(sapply(huhn, as.numeric, simplify = FALSE))
          if (!is.null(sd.max)) {
            huhn <- huhn[huhn$sd < sd.max, ]
          }
          x <- try(huhn[huhn$rt > 0, ], silent = TRUE)
          if (inherits(x, "try-error")) {
            NA
          } else {
            x
          }
        }
        pks <- fit_peaks(
          chrom_list[[sample]],
          lambda = lambda,
          fit = fit,
          max.iter = max.iter,
          sd.max = sd.max,
          estimate_purity = estimate_purity,
          noise_threshold = noise_threshold,
          ...
        )
        pks <- cbind(sample = names(chrom_list)[sample], lambda, pks)
        pks <- remove_bad_peaks(pks)
        pks <- convert_indices_to_times(
          pks,
          chrom_list = chrom_list,
          idx = sample,
          tfac = tfac
        )
        pks
      })
    )
    names(ptable) <- lambdas
    if (collapse) {
      ptable <- do.call(rbind, ptable)
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
    class = "peak_list"
  )
}
