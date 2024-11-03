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
    df <- tidytable::bind_rows(list, .id = "id")
    df <- df |>
      tidytable::group_by(id) |>
      tidytable::mutate(time = time + shift)

    if (normalize_intensity) {
      df <- df |>
        tidytable::mutate(intensity = intensity / max(intensity))
    }
    if (normalize_time) {
      df <- df |>
        tidytable::mutate(time_2 = max(time)) |>
        tidytable::mutate(time = time / time_2) |>
        tidytable::ungroup()
    }
    return(df)
  }
