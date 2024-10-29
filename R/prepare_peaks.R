#' Prepare peaks
#'
#' @param x X
#'
#' @return Prepared peaks
#'
#' @examples NULL
prepare_peaks <- function(x) {
  MSnbase::Chromatogram(
    intensity = x$intensity,
    rtime = (x$time - min(x$time)) /
      (max(x$time) - min(x$time))
  ) ## see https://github.com/sneumann/xcms/issues/593
}
