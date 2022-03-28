#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
extract_ms_peak <- function(x) {
  feature <- seq_along(x)
  y <- mclapply(X = feature, function(z) {
    chrom <- MSnbase::Chromatogram(
      intensity = x[[z]]$intensity,
      rtime = x[[z]]$rtime
    )
    return(chrom)
  })
  return(y)
}
