#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
extract_ms_progress <- function(xs) {
  # p <- progressor(along = xs)
  pbapply::pblapply(
    # lapply(
    X = xs,
    FUN = function(x) {
      # p()
      lapply(
        X = seq_along(
          1:nrow(peaks_prelist$list_df_features_with_peaks_long[[x]])
        ),
        FUN = function(z) {
          ms_chr <-
            ## avoid confusion with MSnbase or xcms
            mzR::chromatogram(
              object = dda_data,
              rt = peaks_prelist$list_rtr[[x]],
              mz = peaks_prelist$list_mzr[[x]][[z]]
            )
          return(ms_chr)
        }
      )
    }
  )
}
