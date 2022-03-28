#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
prepare_peaks <- function(x) {
  MSnbase::Chromatogram(
    intensity = x$intensity,
    rtime = (x$time - min(x$time)) /
      (max(x$time) - min(x$time))
  ) #' see https://github.com/sneumann/xcms/issues/593
}
