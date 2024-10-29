#' Keep best candidates
#'
#' @param df Dataframe
#'
#' @return A dataframe containing the best candidates only
#'
#' @examples NULL
keep_best_candidates <- function(df) {
  best_candidates <- df |>
    dplyr::mutate_all(list(~ gsub(
      pattern = "\\|.*",
      replacement = "",
      x = .x
    ))) |>
    dplyr::mutate(
      best_candidate_1 = gsub(
        pattern = " or .*",
        replacement = "",
        x = candidate_structure_tax_npc_01pat
      ),
      best_candidate_2 = gsub(
        pattern = " or .*",
        replacement = "",
        x = candidate_structure_tax_npc_02sup
      ),
      best_candidate_3 = gsub(
        pattern = " or .*",
        replacement = "",
        x = candidate_structure_tax_npc_03cla
      )
    ) |>
    dplyr::distinct(
      feature_id,
      mz = feature_mz,
      rt = feature_rt,
      smiles_2D = candidate_structure_smiles_no_stereo,
      inchikey_2D = candidate_structure_inchikey_no_stereo,
      molecular_formula = candidate_structure_molecular_formula,
      score_biological,
      score_chemical,
      score_final,
      best_candidate_organism = candidate_structure_organism_occurrence_closest,
      consensus_1 = feature_pred_tax_npc_01pat_val,
      consensus_2 = feature_pred_tax_npc_02sup_val,
      consensus_3 = feature_pred_tax_npc_03cla_val,
      consistency_1 = feature_pred_tax_npc_01pat_score,
      consistency_2 = feature_pred_tax_npc_02sup_score,
      consistency_3 = feature_pred_tax_npc_03cla_score,
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
