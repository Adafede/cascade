#' Add peak metadata
#'
#' @param df Dataframe
#' @param min_score Min score
#'
#' @return A dataframe with peak metadata
#'
#' @examples NULL
add_peak_metadata <- function(df, min_score = 0.8) {
  df_meta <- df |>
    tidytable::arrange(tidytable::desc(feature_area)) |>
    tidytable::distinct(
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
      ## add consensus
      sample,
      keep
    ) |>
    tidytable::group_by(sample, peak_id) |>
    tidytable::distinct(feature_id, .keep_all = TRUE) |>
    tidytable::add_count(name = "featuresPerPeak") |>
    tidytable::distinct(molecular_formula, smiles_2D, inchikey_2D, .keep_all = TRUE) |>
    tidytable::add_count(name = "molecularFormulasPerPeak") |>
    tidytable::distinct(smiles_2D, inchikey_2D, .keep_all = TRUE) |>
    tidytable::add_count(name = "structuresPerPeak") |>
    tidytable::distinct(best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      .keep_all = TRUE
    ) |>
    tidytable::add_count(name = "chemicalClassesPerPeak") |>
    tidytable::distinct(best_candidate_1, best_candidate_2, .keep_all = TRUE) |>
    tidytable::add_count(name = "chemicalSuperclassesPerPeak") |>
    tidytable::distinct(best_candidate_1, .keep_all = TRUE) |>
    tidytable::add_count(name = "chemicalPathwaysPerPeak") |>
    tidytable::ungroup() |>
    tidytable::distinct(
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
    tidytable::arrange(tidytable::desc(feature_area)) |>
    tidytable::filter(comparison_score >= min_score) |>
    tidytable::distinct(
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
      ## add consensus
      sample,
      keep
    ) |>
    tidytable::group_by(sample, peak_id) |>
    tidytable::distinct(feature_id, .keep_all = TRUE) |>
    tidytable::add_count(name = "correlatedFeaturesPerPeak") |>
    tidytable::distinct(molecular_formula, smiles_2D, inchikey_2D, .keep_all = TRUE) |>
    tidytable::add_count(name = "correlatedMolecularFormulasPerPeak") |>
    tidytable::distinct(smiles_2D, inchikey_2D, .keep_all = TRUE) |>
    tidytable::add_count(name = "correlatedStructuresPerPeak") |>
    tidytable::distinct(best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      .keep_all = TRUE
    ) |>
    tidytable::add_count(name = "correlatedChemicalClassesPerPeak") |>
    tidytable::distinct(best_candidate_1, best_candidate_2, .keep_all = TRUE) |>
    tidytable::add_count(name = "correlatedChemicalSuperclassesPerPeak") |>
    tidytable::distinct(best_candidate_1, .keep_all = TRUE) |>
    tidytable::add_count(name = "correlatedChemicalPathwaysPerPeak") |>
    tidytable::ungroup() |>
    tidytable::distinct(
      sample,
      peak_id,
      correlatedFeaturesPerPeak,
      correlatedMolecularFormulasPerPeak,
      correlatedStructuresPerPeak,
      correlatedChemicalClassesPerPeak,
      correlatedChemicalSuperclassesPerPeak,
      correlatedChemicalPathwaysPerPeak
    )

  df_meta_taxo <- df |>
    tidytable::arrange(tidytable::desc(feature_area)) |>
    tidytable::filter(keep == "Y") |>
    tidytable::filter(comparison_score >= min_score) |>
    tidytable::distinct(
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
      ## add consensus
      sample,
      keep
    ) |>
    tidytable::group_by(sample, peak_id) |>
    tidytable::distinct(feature_id, .keep_all = TRUE) |>
    tidytable::add_count(name = "finalFeaturesPerPeak") |>
    tidytable::distinct(molecular_formula, smiles_2D, inchikey_2D, .keep_all = TRUE) |>
    tidytable::add_count(name = "finalMolecularFormulasPerPeak") |>
    tidytable::distinct(smiles_2D, inchikey_2D, .keep_all = TRUE) |>
    tidytable::add_count(name = "finalStructuresPerPeak") |>
    tidytable::distinct(best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      .keep_all = TRUE
    ) |>
    tidytable::add_count(name = "finalChemicalClassesPerPeak") |>
    tidytable::distinct(best_candidate_1, best_candidate_2, .keep_all = TRUE) |>
    tidytable::add_count(name = "finalChemicalSuperclassesPerPeak") |>
    tidytable::distinct(best_candidate_1, .keep_all = TRUE) |>
    tidytable::add_count(name = "finalChemicalPathwaysPerPeak") |>
    tidytable::ungroup() |>
    tidytable::distinct(
      sample,
      peak_id,
      finalFeaturesPerPeak,
      finalMolecularFormulasPerPeak,
      finalStructuresPerPeak,
      finalChemicalClassesPerPeak,
      finalChemicalSuperclassesPerPeak,
      finalChemicalPathwaysPerPeak
    )

  my_bins <- function(x) {
    tidytable::case_when(
      x == 0 ~ "0",
      x == 1 ~ "01",
      x > 1 & x <= 5 ~ "02-05",
      x > 5 & x <= 10 ~ "06-10",
      x > 10 ~ "10+"
    )
  }

  df_meta_full <- df_meta |>
    tidytable::full_join(df_meta_cor) |>
    tidytable::full_join(df_meta_taxo) |>
    tidytable::fill(sample, .direction = "downup")

  df_meta_full[is.na(df_meta_full)] <- 0

  df_meta_full <-
    cbind(df_meta_full, setNames(df_meta_full[colnames(df_meta_full)[grepl(pattern = "PerPeak", x = colnames(df_meta_full))]], paste0(colnames(df_meta_full)[grepl(pattern = "PerPeak", x = colnames(df_meta_full))], "_old")))

  df_meta_full <- df_meta_full |>
    dplyr::mutate_at(.vars = colnames(df_meta_full)[grepl(pattern = "PerPeak$", x = colnames(df_meta_full))], .funs = my_bins)

  return(df_meta_full)
}
