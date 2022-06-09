#' Title
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
baseline_chromatogram <- function(df) {
  df2 <- df
  
  intensity <- df$intensity
  
  intensity[is.na(intensity)] <- 0
  
  intensity_baseline <- baseline(spectra = t(intensity),
                                 method = "peakDetection")
  
  intensity_new <- t(intensity_baseline@corrected) |>
    data.table::data.table()
  
  df2$intensity <- intensity_new$V1
  
  return(df2)
}