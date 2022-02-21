#' Title
#'
#' @param ms_peak
#'
#' @return
#' @export
#'
#' @examples
compare_peaks <-
  function(ms_peak) {
    if (length(ms_peaks[[ms_peak]]) > 1) {
      MSnbase::compareChromatograms(cad_peak,
                                    ms_peaks[[ms_peak]],
                                    ALIGNFUNARGS = list(method = "approx"))
    } else {
      0
    }
  }
