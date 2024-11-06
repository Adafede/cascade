#' Improve signals progress
#'
#' @include improve_signal.R
#'
#' @param xs XS
#' @param fourier_components Fourier components
#' @param frequency Frequency
#' @param resample Resample
#' @param time_min Time min
#' @param time_max Time max
#'
#' @return A list of data frames with improved signals
#'
#' @examples NULL
improve_signals_progress <- function(xs,
                                     fourier_components = 0.01,
                                     frequency = 2,
                                     resample = 1,
                                     time_min = 0,
                                     time_max = Inf) {
  p <- progressr::progressor(along = xs)
  xs |>
    furrr::future_map(
      .f = function(x,
                    fourier_components,
                    frequency,
                    resample,
                    time_min,
                    time_max) {
        p()
        improve_signal(
          df = x |>
            tidytable::select(time, intensity),
          fourier_components = fourier_components,
          frequency = frequency,
          resample = resample,
          time_min = time_min,
          time_max = time_max
        )
      },
      fourier_components = fourier_components,
      frequency = frequency,
      resample = resample,
      time_min = time_min,
      time_max = time_max
    )
}
