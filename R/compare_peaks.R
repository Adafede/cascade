#' Compare peaks
#'
#' @param x X
#'
#' @return A comparison score
#'
#' @examples Null
compare_peaks <- function(x) {
  if (length(list_ms_peaks[[x]]) != 0) {
    future.apply::future_lapply(
      X = seq_along(list_ms_peaks[[x]]),
      FUN = function(z) {
        if (length(list_ms_peaks[[x]][[z]]) > 1) {
          score <- MSnbase::compareChromatograms(
            x = switch(detector,
              "bpi" = peaks_prelist_bpi$list_chromato_peaks,
              "cad" = peaks_prelist_cad$list_chromato_peaks,
              "pda" = peaks_prelist_pda$list_chromato_peaks
            )[[x]],
            y = list_ms_peaks[[x]][[z]],
            method = "closest"
          )
        } else {
          score <- 0
        }
        return(score)
      }
    )
  } else {
    rep(list(0), length(x))
  }
}
