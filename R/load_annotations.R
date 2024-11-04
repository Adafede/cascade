#' Load annotations
#'
#' @param file File
#' @param show_example Show example? Default to FALSE
#' @param mode Mode
#'
#' @return A table of annotations
#'
#' @examples NULL
load_annotations <- function(file = NULL,
                             show_example = FALSE,
                             mode = "pos") {
  if (show_example) {
    # annotation_table |>
    #   saveRDS(file = "inst/extdata/annotations.rds")
    readRDS(system.file("extdata", "annotations.rds", package = "cascade"))
  } else {
    file |>
      tidytable::fread() |>
      tidytable::mutate(tidytable::across(
        tidytable::contains("candidate_structure_tax"),
        .fns = function(x) {
          gsub(
            pattern = "\\$",
            replacement = "or",
            x = x
          )
        }
      )) |>
      tidytable::mutate(tidytable::across(
        tidytable::contains("candidate_structure_tax"),
        .fns = function(x) {
          gsub(
            pattern = "u00A7",
            replacement = "$",
            x = x
          )
        }
      )) |>
      tidytable::mutate(mode = mode)
  }
}
