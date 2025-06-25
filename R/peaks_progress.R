#' Peaks progress
#'
#' @include get_peaks.R
#'
#' @param df_xy Df X Y
#'
#' @return A list of peaks
#'
#' @examples NULL
peaks_progress <- function(df_xy) {
  list(
    "666" = df_xy |>
      tidytable::filter(rtime >= 0) |>
      tidytable::select(rtime, intensity) |>
      tidytable::rename(`666` = intensity) |>
      tibble::column_to_rownames("rtime") |>
      as.matrix()
  ) |>
    get_peaks(
      lambdas = c("666"),
      # fit = c("egh", "gaussian", "raw"),
      # fit = c("egh"),
      # fit = c("gaussian"),
      # fit = c("raw"),
      sd.max = 50,
      max.iter = 1000,
      noise_threshold = 0.001
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
