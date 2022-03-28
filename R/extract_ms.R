#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
extract_ms <- function(x) {
  feature <- seq_along(1:nrow(list_df_peaks_split_flat[[x]]))
  y <- mclapply(X = feature, function(z) {
    ms_chr <- list_dda_split[[x]] |>
      MSnbase::filterRt(c(list_rtr[[x]][[z]][1], list_rtr[[x]][[z]][2])) |>
      MSnbase::filterMz(c(list_mzr[[x]][[z]][1], list_mzr[[x]][[z]][2])) |>
      MSnbase::chromatogram(
        rt = list_rtr[[x]][[z]],
        mz = list_mzr[[x]][[z]]
      )
    return(ms_chr)
  })
  return(y)
}
