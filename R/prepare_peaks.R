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
    rtime = (x$rtime - min(x$rtime)) /
      (max(x$rtime) - min(x$rtime))
  ) ## see https://github.com/sneumann/xcms/issues/593
}
