#' Prepare rt
#'
#' @param x X
#'
#' @return Prepared RTs
#'
#' @examples NULL
prepare_rt <- function(x) {
  rtr <- x[1, ] |>
    tidytable::mutate(
      rtmin = (rt_min + CAD_SHIFT) * 60,
      rtmax = (rt_max + CAD_SHIFT) * 60
    ) |>
    tidytable::select(rtmin, rtmax) |>
    as.matrix()
  return(rtr)
}
