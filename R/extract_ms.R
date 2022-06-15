#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
extract_ms <- function(x) {
  feature <- seq_along(1:nrow(switch(detector,
    "bpi" = peaks_prelist_bpi$list_df_features_with_peaks_long,
    "cad" = peaks_prelist_cad$list_df_features_with_peaks_long,
    "pda" = peaks_prelist_pda$list_df_features_with_peaks_long
  )[[x]]))
  y <- mclapply(X = feature, function(z) {
    ms_chr <- switch(detector,
      "bpi" = peaks_prelist_bpi$list_dda_with_peak,
      "cad" = peaks_prelist_cad$list_dda_with_peak,
      "pda" = peaks_prelist_pda$list_dda_with_peak
    )[[x]] |>
      MSnbase::filterRt(c(switch(detector,
        "bpi" = peaks_prelist_bpi$list_rtr,
        "cad" = peaks_prelist_cad$list_rtr,
        "pda" = peaks_prelist_pda$list_rtr
      )[[x]][[z]][1], switch(detector,
        "bpi" = peaks_prelist_bpi$list_rtr,
        "cad" = peaks_prelist_cad$list_rtr,
        "pda" = peaks_prelist_pda$list_rtr
      )[[x]][[z]][2])) |>
      MSnbase::filterMz(c(switch(detector,
        "bpi" = peaks_prelist_bpi$list_mzr,
        "cad" = peaks_prelist_cad$list_mzr,
        "pda" = peaks_prelist_pda$list_mzr
      )[[x]][[z]][1], switch(detector,
        "bpi" = peaks_prelist_bpi$list_mzr,
        "cad" = peaks_prelist_cad$list_mzr,
        "pda" = peaks_prelist_pda$list_mzr
      )[[x]][[z]][2])) |>
      MSnbase::chromatogram(
        rt = switch(detector,
          "bpi" = peaks_prelist_bpi$list_rtr,
          "cad" = peaks_prelist_cad$list_rtr,
          "pda" = peaks_prelist_pda$list_rtr
        )[[x]][[z]],
        mz = switch(detector,
          "bpi" = peaks_prelist_bpi$list_mzr,
          "cad" = peaks_prelist_cad$list_mzr,
          "pda" = peaks_prelist_pda$list_mzr
        )[[x]][[z]]
      )
    return(ms_chr)
  })
  return(y)
}
