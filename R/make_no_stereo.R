#' Make no stereo
#'
#' @param df Dataframe
#'
#' @return A dataframe with no stereo structures
#'
#' @examples NULL
make_no_stereo <-
  function(df) {
    message("Keeping no stereo structures only")
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
