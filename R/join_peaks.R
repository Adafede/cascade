#' Join peaks
#'
#' @param chromatograms Chromatograms
#' @param peaks Peaks
#' @param min_area Min area
#'
#' @return A dataframe with joined peaks
#'
#' @examples NULL
join_peaks <- function(chromatograms, peaks, min_area) {
  cat("setting joining keys \n")
  data.table::setkey(peaks, rt_min, rt_max)
  data.table::setkey(chromatograms, rt_1, rt_2)

  cat("joining within given rt tolerance \n")
  df <- data.table::foverlaps(peaks, chromatograms) |>
    dplyr::filter(id == i.id) |>
    dplyr::group_by(peak_id, id) |>
    dplyr::mutate(integral = sum(intensity)) |>
    dplyr::ungroup() |>
    dplyr::distinct(peak_id, id, peak_max, rt_apex, rt_min, rt_max, integral) |>
    dplyr::group_by(id) |>
    dplyr::filter(integral / sum(integral) >= min_area) |>
    data.table::data.table()
  return(df)
}
