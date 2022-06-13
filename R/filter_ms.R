#' Title
#'
#' @param x
#' @param shift
#'
#' @return
#' @export
#'
#' @examples
filter_ms <- function(x, shift) {
  MSnbase::filterFile(
    dda_data,
    dda_data@phenoData@data$sampleNames[grepl(
      pattern = unique(x$id),
      x = dda_data@phenoData@data$sampleNames
    )]
  ) |>
    MSnbase::filterRt(rt = c(
      (min(x$rt_min) + shift - RT_TOL) * 60,
      (max(x$rt_max) + shift + RT_TOL) * 60
    )) |>
    MSnbase::filterMz(mz = c(
      min(x$mz_min),
      max(x$mz_max)
    ))
}
