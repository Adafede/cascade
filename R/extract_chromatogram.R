#' Extract chromatogram
#'
#' @include change_intensity_name.R
#'
#' @param list List
#' @param type Type
#'
#' @return An extracted chromatogram
#'
#' @examples NULL
extract_chromatogram <- function(list, type) {
  stopifnot("type must be one of 'bpi, 'cad',or 'pda'" = type %in% c("bpi", "cad", "pda"))
  return(list[switch(type,
    "bpi" = c(TRUE, FALSE, FALSE),
    "cad" = c(FALSE, FALSE, TRUE),
    "pda" = c(FALSE, TRUE, FALSE)
  )] |>
    purrr::pluck(1) |>
    change_intensity_name(name = switch(type,
      "bpi" = "BasePeak_0",
      "cad" = "UV.1_CAD_1_0",
      "pda" = "PDA.1_TotalAbsorbance_0"
    )))
}
