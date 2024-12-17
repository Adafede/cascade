#' Prepare mz
#'
#' @param x X
#'
#' @return A list of prepared mz's
#'
#' @examples NULL
prepare_mz <- function(x) {
  seq_along(seq_len(nrow(x))) |>
    purrr::map(
      .f = function(z) {
        x[z, ] |>
          tidytable::select(mzmin = mz_min, mzmax = mz_max) |>
          as.matrix()
      }
    )
}
