#' Hierarchies grouped progress
#'
#' @param xs XS
#'
#' @return A list of grouped hierarchies
#'
#' @examples NULL
hierarchies_grouped_progress <- function(xs) {
  xs |>
    purrr::map(
      .progress = TRUE,
      .f = function(x) {
        if (nrow(x) != 0) {
          ## dirty workaround
          if (
            nrow(
              x |>
                tidytable::filter(!grepl(pattern = "var. ", x = taxaLabels)) |>
                tidytable::distinct(taxa)
            ) !=
              1
          ) {
            prepare_hierarchy(
              dataframe = x |>
                tidytable::mutate(
                  best_candidate_1 = chemical_pathway,
                  best_candidate_2 = chemical_superclass,
                  best_candidate_3 = chemical_class,
                  organism = taxaLabels
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
          }
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
