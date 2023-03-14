#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
prepare_rt <- function(x) {
  rtr <- x[1, ] |>
    mutate(
      rtmin = (rt_min + CAD_SHIFT) * 60,
      rtmax = (rt_max + CAD_SHIFT) * 60
    ) |>
    select(rtmin, rtmax) |>
    as.matrix()
  return(rtr)
}
