#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
transform_baseline <- function(x) {
  for (i in seq_along(seq_len(length(x)))) {
    intensity <- x[[i]]$intensity

    intensity[is.na(intensity)] <- 0

    intensity_baseline <- baseline(
      spectra = t(intensity),
      method = "peakDetection"
    )

    intensity_new <- t(intensity_baseline@corrected) |>
      data.table::data.table()

    newlist <- x
    newlist[[i]]$intensity <- intensity_new$V1
    return(newlist)
  }
}
