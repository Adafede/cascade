#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
hierarchies_grouped_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      if (nrow(x) != 0) {
        ## dirty workaround
        if (nrow(x |>
          dplyr::filter(!grepl(pattern = "var. ", x = taxaLabels)) |>
          dplyr::distinct(taxa)) != 1) {
          prepare_hierarchy(
            dataframe = x |>
              dplyr::mutate(
                best_candidate_1 = chemical_pathway,
                best_candidate_2 = chemical_superclass,
                best_candidate_3 = chemical_class,
                organism = taxaLabels
              ) |>
              dplyr::mutate(sample = organism, species = organism) |>
              dplyr::select(-taxa, -taxaLabels, -references, -referencesLabels) |>
              dplyr::distinct(),
            type = "literature"
          )
        }
      } else {
        data.frame() |>
          dplyr::mutate(
            parents = NA,
            ids = NA,
            labels = NA,
            sample = NA,
            species = NA,
            values = 0
          ) |>
          dplyr::mutate_if(is.logical, as.character)
      }
    }
  )
}
