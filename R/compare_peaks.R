#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
compare_peaks <- function(x) {
  if (length(list_ms_peaks[[x]]) != 0) {
    feature <- seq_along(list_ms_peaks[[x]])
    y <- mclapply(
      X = feature,
      FUN = function(z) {
        if (length(list_ms_peaks[[x]][[z]]) > 1) {
          score <- MSnbase::compareChromatograms(list_chromato_peaks[[x]],
            list_ms_peaks[[x]][[z]],
            ALIGNFUNARGS = list(method = "approx")
          )
        } else {
          score <- 0
        }
        return(score)
      }
    )
    return(y)
  } else {
    rep(list(0), length(x))
  }
}
