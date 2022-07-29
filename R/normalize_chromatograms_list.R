#' Title
#'
#' @param list
#' @param shift
#' @param time
#' @param intensity
#'
#' @return
#' @export
#'
#' @examples
normalize_chromatograms_list <-
  function(list,
           shift = 0,
           time = FALSE,
           intensity = TRUE) {
    df <- dplyr::bind_rows(list, .id = "id")
    if (intensity == TRUE) {
      df <- df |>
        dplyr::mutate(intensity = intensity / max(intensity))
    }
    if (time == TRUE) {
      df <- df |>
        dplyr::group_by(id) |>
        dplyr::mutate(time = time + shift) |>
        dplyr::mutate(time_2 = max(time)) |>
        dplyr::mutate(time = time / time_2) |>
        dplyr::ungroup()
    }
    return(df)
  }
