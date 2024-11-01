#' Extract MS progress
#'
#' @param xs XS
#'
#' @return A list of extracted MS peaks
#'
#' @examples NULL
extract_ms_progress <- function(xs) {
  # p <- progressor(along = xs)
  pbapply::pblapply(
    # lapply(
    X = xs,
    FUN = function(x) {
      # p()
      message("CAD Peak: ", x)
      lapply(
        X = seq_along(
          seq_len(nrow(peaks_prelist$list_df_features_with_peaks_long[[x]]))
        ),
        FUN = function(z) {
          # message("CAD Peak: ", x, ", MS feature: ", z)
          tryCatch(
            expr = {
              dda_data |>
                MSnbase::chromatogram(rt = peaks_prelist$list_rtr[[x]], mz = peaks_prelist$list_mzr[[x]][[z]])
            },
            error = function(e) {
              message("Going too fast...1 sec")
              Sys.sleep(1)
              tryCatch(
                expr = {
                  dda_data |>
                    MSnbase::chromatogram(rt = peaks_prelist$list_rtr[[x]], mz = peaks_prelist$list_mzr[[x]][[z]])
                },
                error = function(e) {
                  message("Going way too fast...2 secs")
                  Sys.sleep(2)
                  return(
                    dda_data |>
                      MSnbase::chromatogram(rt = peaks_prelist$list_rtr[[x]], mz = peaks_prelist$list_mzr[[x]][[z]])
                  )
                }
              )
            }
          )
        }
      )
    }
  )
}
