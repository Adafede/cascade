#' Peaks progress
#'
#' @param xs XS
#'
#' @return A list of peaks
#'
#' @examples NULL
peaks_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      # plot(cads_baselined$intensity, type = "l", col = "navy")
      # grid()
      found <-
        pracma::findpeaks(x$intensity,
          threshold = 0.001,
          sortstr = TRUE
        )
      # points(x[, 2], x[, 1], pch = 20, col = "maroon") ## End(Not run)

      peaks <- data.frame(found) |>
        dplyr::mutate(
          peak_id = dplyr::row_number(),
          peak_max = X1,
          rt_apex = x$time[X2],
          rt_min = x$time[X3],
          rt_max = x$time[X4]
        )

      data_cascade <- list()
      data_cascade[[1]] <-
        x |>
        dplyr::filter(time >= 0) |>
        tibble::column_to_rownames("time") |>
        dplyr::select(`666` = intensity) |>
        as.matrix()
      names(data_cascade)[[1]] <- "666"

      data <-
        chromatographR::preprocess(
          X = data_cascade,
          dim1 = seq(
            from = as.numeric(rownames(data_cascade[[1]])[[1]]),
            to = as.numeric(tail(rownames(data_cascade[[1]]), 1)),
            by = (abs(
              as.numeric(rownames(data_cascade[[1]])[[2]]) -
                as.numeric(rownames(data_cascade[[1]])[[1]])
            ))
          ),
          dim2 = "666",
          remove.time.baseline = FALSE,
          spec.smooth = FALSE
        )

      pks <- chromatographR::get_peaks(
        chrom_list = data,
        lambdas = c("666"),
        # fit = c("egh", "gaussian", "raw"),
        # fit = c("egh"),
        # fit = c("gaussian"),
        # fit = c("raw"),
        sd.max = 50,
        max.iter = 1000,
        noise_threshold = 0.001,
        show_progress = TRUE
      )

      peaks <- pks$`666`$`666` |>
        dplyr::mutate(
          peak_id = dplyr::row_number(),
          peak_max = height,
          rt_apex = rt,
          rt_min = start,
          rt_max = end
        )

      return(peaks)
    }
  )
}
