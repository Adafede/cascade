#' Title
#'
#' @param list 
#' @param shift 
#'
#' @return
#' @export
#'
#' @examples
normalize_chromatograms_list <- function(list, shift = 0) {
  df <- dplyr::bind_rows(list, .id = "id") |>
    dplyr::mutate(intensity = intensity / max(intensity),
                  time = time + shift) |>
    dplyr::group_by(id) |>
    dplyr::mutate(time_2 = max(time)) |>
    dplyr::mutate(time = time / time_2)
  return(df)
}
