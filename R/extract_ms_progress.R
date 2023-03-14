#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
extract_ms_progress <- function(xs) {
  p <- progressor(along = xs)
  future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      lapply(
        X = seq_along(
          1:nrow(peaks_prelist$list_df_features_with_peaks_long[[x]])
        ),
        FUN = function(z) {
          ms_chr <- dda_data |>
            chromatogram(
              rt = peaks_prelist$list_rtr[[x]],
              mz = peaks_prelist$list_mzr[[x]][[z]]
            )
          return(ms_chr)
        }
      )
    }
  )
}
