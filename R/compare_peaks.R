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
    MSnbase::compareChromatograms(cad_peak,
      ms_peak,
      ALIGNFUNARGS = list(method = "approx")
    )
  }
