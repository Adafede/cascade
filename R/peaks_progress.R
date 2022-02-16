#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
peaks_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      # plot(cads_baselined$intensity, type = "l", col = "navy")
      # grid()
      found <-
        pracma::findpeaks(
          x$intensity,
          npeaks = 2000,
          threshold = 0.01,
          sortstr = TRUE
        )
      # points(x[, 2], x[, 1], pch = 20, col = "maroon") ## End(Not run)
      
      peaks <- data.frame(found) |>
        dplyr::mutate(
          peak_id = dplyr::row_number(),
          peak_max = X1,
          rt_apex = chromatograms_cad_baselined[[i]]$time[X2],
          rt_min = chromatograms_cad_baselined[[i]]$time[X3],
          rt_max = chromatograms_cad_baselined[[i]]$time[X4]
        ) |>
        data.table::data.table()
      
      return(peaks)
    }
  )
}
