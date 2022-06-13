#' Title
#'
#' @param x
#' @param detector
#'
#' @return
#' @export
#'
#' @examples
extract_ms <- function(x, detector) {
  list_df_peaks_split_flat <- switch(detector,
    "cad" = list_df_peaks_cad_split_flat,
    "pda" = list_df_peaks_pda_split_flat
  )
  list_dda_split <- switch(detector,
    "cad" = list_dda_split_cad,
    "pda" = list_dda_split_pda
  )
  list_rtr <- switch(detector,
    "cad" = list_rtr_cad,
    "pda" = list_rtr_pda
  )
  list_mzr <- switch(detector,
    "cad" = list_mzr_cad,
    "pda" = list_mzr_pda
  )

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
