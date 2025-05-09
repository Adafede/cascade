#' Extract MS peak
#'
#' @param x X
#'
#' @return A peak
#'
#' @examples NULL
extract_ms_peak <- function(x) {
  seq_along(x) |>
    purrr::map(
      .f = function(z) {
        MSnbase::Chromatogram(
          intensity = x[[z]]$intensity,
          rtime = x[[z]]$rtime
        )
      }
    )
}
