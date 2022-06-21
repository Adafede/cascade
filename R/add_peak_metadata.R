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
    dplyr::arrange(desc(intensity)) |>
    dplyr::distinct(
      id,
      peak_id,
      integral,
      feature_id,
      rt,
      mz,
      smiles_2D,
      inchikey_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      score_biological,
      score_chemical,
      score_final,
      #' add consensus
      sample,
      species
    ) |>
    dplyr::group_by(id, peak_id) |>
    dplyr::distinct(feature_id,
      # smiles_2D,
      # inchikey_2D,
      # best_candidate_1,
      # best_candidate_2,
      # best_candidate_3,
      .keep_all = TRUE
    ) |>
    dplyr::add_count(name = "featuresPerPeak") |>
    dplyr::distinct(smiles_2D,
      inchikey_2D,
      # best_candidate_1,
      # best_candidate_2,
      # best_candidate_3,
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
      id,
      peak_id,
      featuresPerPeak,
      structuresPerPeak,
      chemicalClassesPerPeak,
      chemicalSuperclassesPerPeak,
      chemicalPathwaysPerPeak
    )

  return(df_meta)
}
