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
      tidytable::group_by(!!as.name(names(df)[!grepl(pattern = "structure", x = names(df))])) |>
      tidytable::distinct(tidytable::any_of(
        c(
          "structure_smiles_2D",
          "reference_wikidata",
          "reference_doi"
        )
      ), .keep_all = TRUE) |>
      tidytable::ungroup() |>
      tidytable::select(-structure_smiles_2D)
  }
