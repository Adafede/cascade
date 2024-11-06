#' Peaks progress
#'
#' @param xs XS
#'
#' @return A list of peaks
#'
#' @examples NULL
peaks_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  xs |>
    furrr::future_map(
      .f = function(x) {
        p()
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
            spec.smooth = FALSE,
            interpolate_rows = FALSE,
            interpolate_cols = FALSE
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
