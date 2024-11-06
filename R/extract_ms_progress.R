#' Extract MS progress
#'
#' @param xs XS
#' @param ms_data MS Data
#' @param peaks_prelist Peaks prelist
#'
#' @return A list of extracted MS peaks
#'
#' @examples NULL
extract_ms_progress <- function(xs, ms_data, peaks_prelist) {
  p <- progressr::progressor(along = xs)
  xs |>
    furrr::future_map(
      .f = function(x, ms_data, peaks_prelist) {
        p()
        message("CAD Peak: ", x)
        lapply(
          X = seq_along(seq_len(
            nrow(peaks_prelist$list_df_features_with_peaks_long[[x]])
          )),
          FUN = function(z, ms_data, peaks_prelist) {
            # message("CAD Peak: ", x, ", MS feature: ", z)
            tryCatch(
              expr = {
                ms_data |>
                  MSnbase::chromatogram(rt = peaks_prelist$list_rtr[[x]], mz = peaks_prelist$list_mzr[[x]][[z]])
              },
              error = function(e) {
                message("Going too fast...1 sec")
                Sys.sleep(1)
                tryCatch(
                  expr = {
                    ms_data |>
                      MSnbase::chromatogram(rt = peaks_prelist$list_rtr[[x]], mz = peaks_prelist$list_mzr[[x]][[z]])
                  },
                  error = function(e) {
                    message("Going way too fast...2 secs")
                    Sys.sleep(2)
                    return(
                      ms_data |>
                        MSnbase::chromatogram(rt = peaks_prelist$list_rtr[[x]], mz = peaks_prelist$list_mzr[[x]][[z]])
                    )
                  }
                )
              }
            )
          },
          ms_data = ms_data,
          peaks_prelist = peaks_prelist
        )
      },
      ms_data = ms_data,
      peaks_prelist = peaks_prelist
    )
}
