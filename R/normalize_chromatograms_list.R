#' Normalize chromatograms list
#'
#' @param list List
#' @param shift Shift
#' @param normalize_intensity Normalize time
#' @param normalize_time Normalize intensity
#'
#' @return A dataframe with normalized chromatograms
#'
#' @examples NULL
normalize_chromatograms_list <-
  function(list,
           shift = 0,
           normalize_intensity = TRUE,
           normalize_time = FALSE) {
    df <- dplyr::bind_rows(list, .id = "id")
    df <- df |>
      dplyr::group_by(id) |>
      dplyr::mutate(time = time + shift)

    if (normalize_intensity) {
      df <- df |>
        dplyr::mutate(intensity = intensity / max(intensity))
    }
    if (normalize_time) {
      df <- df |>
        dplyr::mutate(time_2 = max(time)) |>
        dplyr::mutate(time = time / time_2) |>
        dplyr::ungroup()
    }
    return(df)
  }
