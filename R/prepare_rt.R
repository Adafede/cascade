#' Prepare rt
#'
#' @param x X
#'
#' @return Prepared RTs
#'
#' @examples NULL
prepare_rt <- function(x, shift = 0) {
  rtr <- x[1, ] |>
    tidytable::mutate(
      rtmin = (rt_min + shift) * 60,
      rtmax = (rt_max + shift) * 60
    ) |>
    tidytable::select(rtmin, rtmax) |>
    as.matrix()
  return(rtr)
}
