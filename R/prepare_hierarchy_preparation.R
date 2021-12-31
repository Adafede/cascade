require(package = dplyr, quietly = TRUE)
require(package = splitstackshape, quietly = TRUE)

#' Title
#'
#' @param dataframe
#'
#' @return
#'
#' @examples
prepare_hierarchy_preparation <- function(dataframe) {
  ms1_best_candidate <- dataframe |>
    dplyr::mutate_all(list(~ gsub(
      pattern = "\\|.*",
      replacement = "",
      x = .
    ))) |>
    splitstackshape::cSplit("best_candidate", sep = "§") |>
    dplyr::mutate(id = 1) |>
    dplyr::distinct(
      id,
      feature_id,
      smiles_2D,
      inchikey_2D,
      score_biological,
      score_chemical,
      score_final,
      best_candidate_organism,
      # consensus_1 = consensus_pat,
      # consensus_2 = consensus_sup,
      # consensus_3 = consensus_cla,
      # consistency_1 = consistency_pat,
      # consistency_2 = consistency_sup,
      # consistency_3 = consistency_cla,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3
    ) |>
    dplyr::mutate_all(list(~ y_as_na(x = ., y = ""))) |>
    dplyr::mutate(
      best_candidate_1 = if_else(
        condition = is.na(smiles_2D),
        true = "notAnnotated",
        false = best_candidate_1
      ),
      best_candidate_2 = if_else(
        condition = is.na(smiles_2D),
        true = paste(best_candidate_1, "notAnnotated"),
        false = best_candidate_2
      ),
      best_candidate_3 = if_else(
        condition = is.na(smiles_2D),
        true = paste(best_candidate_2, "notAnnotated"),
        false = best_candidate_3
      ),
      best_candidate_1 = if_else(
        condition = !is.na(smiles_2D) &
          is.na(best_candidate_1),
        true = "notClassified",
        false = best_candidate_1
      ),
      best_candidate_2 = if_else(
        condition = !is.na(smiles_2D) &
          is.na(best_candidate_2),
        true = paste(best_candidate_1, "notClassified"),
        false = best_candidate_2
      ),
      best_candidate_3 = if_else(
        condition = !is.na(smiles_2D) &
          is.na(best_candidate_3),
        true = paste(best_candidate_2, "notClassified"),
        false = best_candidate_3
      )
    )

  ms1_multiple <- ms1_best_candidate |>
    dplyr::left_join(top_m) |>
    #' add this step
    dplyr::filter(!is.na(species)) |>
    dplyr::filter(intensity != 0)

  return(ms1_multiple)
}
