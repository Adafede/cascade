#' Peaks progress
#'
#' @include get_peaks.R
#'
#' @param df_xy Df X Y
#' @param sd_max Maximum standard deviation for peak filtering. Default is 50.
#' @param max_iter Maximum iterations for peak fitting. Default is 1000.
#' @param noise_threshold Noise threshold for peak detection. Default is 0.001.
#' @param fit Peak fitting method. One of "egh", "gaussian", or "raw". Default
#'   is "egh".
#'
#' @return A list of peaks
#'
#' @examples NULL
peaks_progress <- function(
  df_xy,
  sd_max = 50,
  max_iter = 1000,
  noise_threshold = 0.001,
  fit = "egh"
) {
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
      fit = fit,
      sd.max = sd_max,
      max.iter = max_iter,
      noise_threshold = noise_threshold
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
