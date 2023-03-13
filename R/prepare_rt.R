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
  y <- mclapply(
    X = feature,
    FUN = function(z) {
      rtr <- x[z, ] |>
        mutate(
          rtmin = (rt_min + CAD_SHIFT) * 60,
          rtmax = (rt_max + CAD_SHIFT) * 60
        ) |>
        select(rtmin, rtmax) |>
        as.matrix()
      return(rtr)
    }
  )
  return(y)
}
