#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
add_peak_metadata <- function(df) {
  df_meta <- df |>
    dplyr::arrange(desc(feature_area)) |>
    dplyr::distinct(
      peak_id,
      feature_id,
      rt,
      mz,
      smiles_2D,
      inchikey_2D,
      molecular_formula,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      score_biological,
      score_chemical,
      score_final,
      #' add consensus
      sample
    ) |>
    dplyr::group_by(sample, peak_id) |>
    dplyr::distinct(feature_id,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "featuresPerPeak") |>
    dplyr::distinct(molecular_formula,
      smiles_2D,
      inchikey_2D,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "molecularFormulasPerPeak") |>
    dplyr::distinct(smiles_2D,
      inchikey_2D,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "structuresPerPeak") |>
    dplyr::distinct(best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "chemicalClassesPerPeak") |>
    dplyr::distinct(best_candidate_1,
      best_candidate_2,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "chemicalSuperclassesPerPeak") |>
    dplyr::distinct(best_candidate_1,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "chemicalPathwaysPerPeak") |>
    dplyr::ungroup() |>
    dplyr::distinct(
      sample,
      peak_id,
      featuresPerPeak,
      molecularFormulasPerPeak,
      structuresPerPeak,
      chemicalClassesPerPeak,
      chemicalSuperclassesPerPeak,
      chemicalPathwaysPerPeak
    )

  df_meta_cor <- df |>
    dplyr::arrange(desc(feature_area)) |>
    dplyr::filter(comparison_score >= params$chromato$peak$similarity$filter) |>
    dplyr::distinct(
      peak_id,
      feature_id,
      rt,
      mz,
      smiles_2D,
      inchikey_2D,
      molecular_formula,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      score_biological,
      score_chemical,
      score_final,
      #' add consensus
      sample
    ) |>
    dplyr::group_by(sample, peak_id) |>
    dplyr::distinct(feature_id,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "correlatedFeaturesPerPeak") |>
    dplyr::distinct(molecular_formula,
      smiles_2D,
      inchikey_2D,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "correlatedMolecularFormulasPerPeak") |>
    dplyr::distinct(smiles_2D,
      inchikey_2D,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "correlatedStructuresPerPeak") |>
    dplyr::distinct(best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "correlatedChemicalClassesPerPeak") |>
    dplyr::distinct(best_candidate_1,
      best_candidate_2,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "correlatedChemicalSuperclassesPerPeak") |>
    dplyr::distinct(best_candidate_1,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "correlatedChemicalPathwaysPerPeak") |>
    dplyr::ungroup() |>
    dplyr::distinct(
      sample,
      peak_id,
      correlatedFeaturesPerPeak,
      correlatedMolecularFormulasPerPeak,
      correlatedStructuresPerPeak,
      correlatedChemicalClassesPerPeak,
      correlatedChemicalSuperclassesPerPeak,
      correlatedChemicalPathwaysPerPeak
    )

  df_meta_full <- df_meta |>
    dplyr::full_join(df_meta_cor) |>
    tidyr::fill(sample, .direction = "downup")

  df_meta_full[is.na(df_meta_full)] <- 0

  return(df_meta_full)
}
