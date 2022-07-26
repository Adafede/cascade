#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
keep_best_candidates <- function(df) {
  best_candidates <- df |>
    dplyr::mutate_all(list(~ gsub(
      pattern = "\\|.*",
      replacement = "",
      x = .x
    ))) |>
    splitstackshape::cSplit("best_candidate", sep = "$") |>
    dplyr::distinct(
      feature_id,
      mz,
      rt,
      smiles_2D,
      inchikey_2D,
      score_biological,
      score_chemical,
      score_final,
      best_candidate_organism,
      consensus_1 = consensus_pat,
      consensus_2 = consensus_sup,
      consensus_3 = consensus_cla,
      consistency_1 = consistency_pat,
      consistency_2 = consistency_sup,
      consistency_3 = consistency_cla,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      mode
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
  return(best_candidates)
}
