#' Make other
#'
#' @param dataframe Dataframe
#' @param value Value
#'
#' @return A dataframe with harmonized "other" subcategories
#'
#' @examples NULL
make_other <- function(dataframe, value = "peak_area") {
  top_4 <- dataframe |>
    tidytable::group_by(best_candidate_1, best_candidate_2) |>
    tidytable::mutate(new = sum(!!as.name(value))) |>
    tidytable::group_by(best_candidate_1) |>
    tidytable::slice_max(new, n = 4, with_ties = FALSE) |>
    tidytable::ungroup() |>
    tidytable::select(
      smiles_2D,
      inchikey_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3
    ) |>
    tidytable::distinct()

  last <- dataframe |>
    tidytable::ungroup() |>
    tidytable::anti_join(
      top_4,
      by = c(
        "best_candidate_1" = "best_candidate_1",
        "best_candidate_2" = "best_candidate_2"
      )
    ) |>
    tidytable::mutate(best_candidate_2 = "Other") |>
    tidytable::select(
      smiles_2D,
      inchikey_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3
    ) |>
    tidytable::distinct()

  new <- tidytable::bind_rows(top_4, last) |>
    tidytable::select(-smiles_2D, -inchikey_2D)

  dataframe |>
    tidytable::select(-best_candidate_2) |>
    tidytable::left_join(new) |>
    tidytable::distinct()
}
