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
        list(
          "666" = x |>
            tidytable::filter(time >= 0) |>
            tidytable::select(time, intensity) |>
            tidytable::rename(`666` = intensity) |>
            tibble::column_to_rownames("time") |>
            as.matrix()
        ) |>
          chromatographR::get_peaks(
            lambdas = c("666"),
            # fit = c("egh", "gaussian", "raw"),
            # fit = c("egh"),
            # fit = c("gaussian"),
            # fit = c("raw"),
            sd.max = 50,
            max.iter = 1000,
            noise_threshold = 0.001,
            show_progress = FALSE,
            cl = "future"
          ) |>
          purrr::pluck("666") |>
          purrr::pluck("666") |>
          tidytable::mutate(
            peak_id = tidytable::row_number(),
            peak_max = height,
            rt_apex = rt,
            rt_min = start,
            rt_max = end
          )
      }
    )
}
