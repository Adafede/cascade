#' Keep best candidates
#'
#' @include y_as_na.R
#'
#' @param df Dataframe
#'
#' @return A dataframe containing the best candidates only
#'
#' @examples NULL
keep_best_candidates <- function(df) {
  df |>
    tidytable::mutate(tidytable::across(
      .cols = tidytable::everything(),
      .fns = ~ gsub(
        pattern = "\\|.*",
        replacement = "",
        x = .x
      )
    )) |>
    tidytable::mutate(
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
    tidytable::distinct(
      feature_id,
      mz = feature_mz,
      rt = feature_rt,
      smiles_2D = candidate_structure_smiles_no_stereo,
      inchikey_2D = candidate_structure_inchikey_connectivity_layer,
      molecular_formula = candidate_structure_molecular_formula,
      candidate_count_similarity_peaks_matched,
      candidate_score_similarity,
      score_initial,
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
    tidytable::mutate(tidytable::across(
      .cols = tidytable::everything(),
      .fns = ~ y_as_na(x = ., y = "")
    )) |>
    tidytable::mutate(
      best_candidate_1 = tidytable::if_else(
        condition = is.na(smiles_2D),
        true = "notAnnotated",
        false = best_candidate_1
      ),
      best_candidate_2 = tidytable::if_else(
        condition = is.na(smiles_2D),
        true = paste(best_candidate_1, "notAnnotated"),
        false = best_candidate_2
      ),
      best_candidate_3 = tidytable::if_else(
        condition = is.na(smiles_2D),
        true = paste(best_candidate_2, "notAnnotated"),
        false = best_candidate_3
      ),
      best_candidate_1 = tidytable::if_else(
        condition = !is.na(smiles_2D) &
          is.na(best_candidate_1),
        true = "notClassified",
        false = best_candidate_1
      ),
      best_candidate_2 = tidytable::if_else(
        condition = !is.na(smiles_2D) &
          is.na(best_candidate_2),
        true = paste(best_candidate_1, "notClassified"),
        false = best_candidate_2
      ),
      best_candidate_3 = tidytable::if_else(
        condition = !is.na(smiles_2D) &
          is.na(best_candidate_3),
        true = paste(best_candidate_2, "notClassified"),
        false = best_candidate_3
      )
    )
}
