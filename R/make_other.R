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
    dplyr::group_by(best_candidate_1, best_candidate_2) |>
    dplyr::mutate(new = sum(!!as.name(value))) |>
    dplyr::group_by(best_candidate_1) |>
    dplyr::slice_max(new, n = 4, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(
      smiles_2D,
      inchikey_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3
    ) |>
    dplyr::distinct()

  last <- dataframe |>
    dplyr::ungroup() |>
    dplyr::anti_join(
      top_4,
      by = c(
        "best_candidate_1" = "best_candidate_1",
        "best_candidate_2" = "best_candidate_2"
      )
    ) |>
    dplyr::mutate(best_candidate_2 = "Other") |>
    dplyr::select(
      smiles_2D,
      inchikey_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3
    ) |>
    dplyr::distinct()

  new <- dplyr::bind_rows(
    top_4,
    last
  ) |>
    dplyr::select(
      -smiles_2D,
      -inchikey_2D
    )

  df_new <- dataframe |>
    dplyr::select(
      # -best_candidate_1, # needed
      -best_candidate_2
    ) |>
    dplyr::left_join(new) |>
    dplyr::distinct()

  return(df_new)
}
