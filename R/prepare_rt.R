#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
prepare_rt <- function(x) {
  feature <- seq_along(1:nrow(x))
  y <- mclapply(X = feature, function(z) {
    rtr <- x[z, ] |>
      dplyr::mutate(
        rtmin = (rt_min + CAD_SHIFT) * 60,
        rtmax = (rt_max + CAD_SHIFT) * 60
      ) |>
      dplyr::select(rtmin, rtmax) |>
      as.matrix()
    return(rtr)
  })
  return(y)
}
