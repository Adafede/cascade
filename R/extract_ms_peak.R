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
  y <- future_lapply(X = feature, function(z) {
    chrom <- Chromatogram(
      intensity = x[[z]]$intensity,
      rtime = x[[z]]$rtime
    )
    return(chrom)
  })
  return(y)
}
