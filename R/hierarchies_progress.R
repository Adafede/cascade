#' Hierarchies Progress
#'
#' @param xs XS
#' @param comparison Comparison
#'
#' @return A list of hierarchies
#'
#' @examples NULL
hierarchies_progress <- function(xs, comparison) {
  if (comparison |> is.null()) {
    comparison <- NA
  }
  xs |>
    purrr::map(
      .progress = TRUE,
      .f = function(x) {
        if (nrow(x) != 0) {
          prepare_hierarchy(
            dataframe = x |>
              tidytable::mutate(
                best_candidate_1 = chemical_pathway,
                best_candidate_2 = chemical_superclass,
                best_candidate_3 = chemical_class
              ) |>
              tidytable::mutate_rowwise(
                organism = tidytable::if_else(
                  condition = i |>
                    grepl(pattern = "\\W+") |>
                    all(),
                  true = taxaLabels,
                  false = taxaLabels |>
                    gsub(pattern = " .*", replacement = ""),
                )
              ) |>
              tidytable::mutate(sample = organism, species = organism) |>
              tidytable::select(
                -taxa,
                -taxaLabels,
                -references,
                -referencesLabels
              ) |>
              tidytable::distinct(),
            type = "literature"
          )
        } else {
          data.frame() |>
            tidytable::mutate(
              parents = NA,
              ids = NA,
              labels = NA,
              sample = NA,
              species = NA,
              values = 0
            ) |>
            tidytable::mutate(tidytable::across(
              .cols = tidytable::where(is.logical),
              .fns = as.character
            ))
        }
      }
    )
}
