#' Make 2D
#'
#' @param df Dataframe
#'
#' @return A dataframe with 2D structures
#'
#' @examples NULL
make_2D <-
  function(df) {
    message("Keeping 2D structures only")
    df |>
      dplyr::group_by(!!as.name(names(df)[!grepl(pattern = "structure", x = names(df))])) |>
      dplyr::distinct(dplyr::across(dplyr::any_of(
        c(
          "structure_smiles_2D",
          "reference_wikidata",
          "reference_doi"
        )
      )), .keep_all = TRUE) |>
      dplyr::ungroup() |>
      dplyr::select(-structure_smiles_2D)
  }
