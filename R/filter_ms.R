#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
filter_ms <- function(x) {
  MSnbase::filterFile(
    dda_data,
    dda_data@phenoData@data$sampleNames[grepl(
      pattern = unique(x$id),
      x = dda_data@phenoData@data$sampleNames
    )]
  ) |>
    MSnbase::filterRt(rt = c(
      (min(x$rt_min) + CAD_SHIFT - RT_TOL) * 60,
      (max(x$rt_max) + CAD_SHIFT + RT_TOL) * 60
    )) |>
    MSnbase::filterMz(mz = c(
      min(x$mz_min),
      max(x$mz_max)
    ))
}
