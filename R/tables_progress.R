#' Tables progress
#'
#' @param xs XS
#' @param structures_classified structures classified
#'
#' @return A list of tables
#'
#' @examples NULL
tables_progress <- function(xs, structures_classified) {
  p <- progressr::progressor(along = xs)
  xs |>
    furrr::future_map(
      .f = function(x, structures_classified) {
        p()
        if (nrow(x != 0)) {
          x |>
            tidytable::left_join(structures_classified) |>
            tidytable::mutate(structureImage = URLencode(structureSmiles)) |>
            tidytable::relocate(structureImage, .after = structure) |>
            tidytable::relocate(structureLabel, .before = structure) |>
            tidytable::select(-references_ids, -structure_id, -structureSmiles) |>
            tidytable::separate_longer_delim(
              c("taxa", "taxaLabels", "references", "referencesLabels"),
              delim = "|"
            ) |>
            tidytable::group_by(structure) |>
            tidytable::fill(c("taxa", "taxaLabels", "references", "referencesLabels"),
              .direction = "downup"
            ) |>
            tidytable::group_by(structureLabel) |>
            tidytable::add_count(sort = TRUE) |>
            tidytable::select(-n) |>
            tidytable::group_by(chemical_class) |>
            tidytable::add_count(sort = TRUE) |>
            tidytable::select(-n) |>
            tidytable::group_by(chemical_superclass) |>
            tidytable::add_count(sort = TRUE) |>
            tidytable::select(-n) |>
            tidytable::group_by(chemical_pathway) |>
            tidytable::add_count(sort = TRUE) |>
            tidytable::select(-n) |>
            tidytable::distinct()
        } else {
          data.frame() |>
            tidytable::mutate(
              structureLabel = NA,
              structure = NA,
              structureImage = NA,
              taxaLabels = NA,
              taxa = NA,
              referenceLabels = NA,
              references = NA,
              chemical_pathway = NA,
              chemical_superclass = NA,
              chemical_class = NA
            )
        }
      },
      structures_classified = structures_classified
    )
}
