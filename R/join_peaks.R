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
  data.table::setkey(peaks, rt_min, rt_max)
  data.table::setkey(chromatograms, rt_1, rt_2)

  data.table::foverlaps(peaks, chromatograms) |>
    # tidytable::filter(id == i.id) |>
    tidytable::group_by(peak_id, id) |>
    tidytable::mutate(integral = sum(intensity)) |>
    tidytable::ungroup() |>
    tidytable::distinct(peak_id, id, peak_max, rt_apex, rt_min, rt_max, integral) |>
    tidytable::group_by(id) |>
    tidytable::filter(integral / sum(integral) >= min_area) |>
    tidytable::data.table()
}
