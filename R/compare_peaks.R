#' Compare peaks
#'
#' @param x X
#' @param list_ms_peaks list_ms_peaks
#' @param peaks_prelist peaks_prelist
#'
#' @return A comparison score
#'
#' @examples NULL
compare_peaks <- function(x, list_ms_peaks, peaks_prelist) {
  if (length(list_ms_peaks[[x]]) != 0) {
    seq_along(list_ms_peaks[[x]]) |>
      purrr::map(
        .f = function(list_ms_peaks, x, z) {
          if (length(list_ms_peaks[[x]][[z]]) > 1) {
            MSnbase::compareChromatograms(peaks_prelist$list_chromato_peaks[[x]],
              y = list_ms_peaks[[x]][[z]],
              method = "closest"
            )
          } else {
            0
          }
        },
        list_ms_peaks = list_ms_peaks,
        x = x
      )
  } else {
    rep(list(0), length(x))
  }
}
