#' Extract MS peak
#'
#' @param x X
#'
#' @return A peak
#'
#' @examples NULL
extract_ms_peak <- function(x) {
  future.apply::future_lapply(X = seq_along(x), function(z) {
    chrom <- MSnbase::Chromatogram(
      intensity = x[[z]]$intensity,
      rtime = x[[z]]$rtime
    )
    return(chrom)
  })
}
