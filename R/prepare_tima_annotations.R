#' Prepare TIMA annotations
#'
#' @export
#'
#' @include keep_best_candidates.R
#' @include load_annotations.R
#' @include make_confident.R
#'
#' @param annotations annotations
#' @param predicted_classes Show predicted classes? Default to FALSE
#' @param min_score_initial Minimal initial score
#' @param min_score_biological Minimal biological score
#' @param min_score_chemical Minimal chemical score
#' @param min_score_final Minimal final score
#' @param min_matched_peaks_absolute Minimal number of matched peaks
#' @param min_matched_peaks_percentage Minimal percentage of matched peaks
#' @param min_peaks Minimal number of peaks in spectrum
#' @param libraries Libraries to consider
#' @param show_example Show example? Default to FALSE
#'
#' @return Prepared tables
#'
#' @examples NULL
prepare_tima_annotations <- function(
  annotations = NULL,
  predicted_classes = FALSE,
  min_score_initial = 0.0,
  min_score_biological = 0.0,
  min_score_chemical = 0.0,
  min_score_final = 0.0,
  min_matched_peaks_absolute = 0L,
  min_matched_peaks_percentage = 0.0,
  min_peaks = 3L,
  libraries = c(
    "gnps",
    "massbank",
    "merlin",
    "ISDB",
    "ISDB - Wikidata",
    "TIMA MS1"
  ),
  show_example = FALSE
) {
  if (annotations |> is.null()) {
    annotations <- c("annotations" = tempfile())
  }
  if (predicted_classes) {
    names(annotations) <- names(annotations) |>
      paste("pred", sep = "_")
  }
  tables_prepared <- annotations |>
    purrr::imap(.f = function(x, name) {
      table <- x |>
        load_annotations(show_example = show_example) |>
        tidytable::mutate(tidytable::across(
          .cols = tidytable::everything(),
          .fns = function(x) {
            x |> gsub(pattern = "\\|.*", replacement = "")
          }
        ))
      if (!predicted_classes) {
        table <- table |>
          tidytable::filter(
            score_initial |> as.numeric() >= min_score_initial
          ) |>
          tidytable::filter(
            candidate_count_similarity_peaks_matched >=
              feature_spectrum_peaks |>
                as.numeric() *
                min_matched_peaks_absolute |
              candidate_count_similarity_peaks_matched |> as.numeric() >=
                min_matched_peaks_absolute
          ) |>
          tidytable::filter(
            feature_spectrum_peaks |> as.numeric() >= min_peaks
          ) |>
          tidytable::filter(candidate_library %in% libraries)
      }
      if (predicted_classes) {
        table <- table |>
          tidytable::mutate(
            candidate_structure_tax_npc_01pat = feature_pred_tax_npc_01pat_val,
            candidate_structure_tax_npc_02sup = feature_pred_tax_npc_02sup_val,
            candidate_structure_tax_npc_03cla = feature_pred_tax_npc_03cla_val
          )
      }
      table <- table |>
        tidytable::filter(
          score_biological |> as.numeric() >= min_score_biological
        ) |>
        tidytable::filter(
          score_chemical |> as.numeric() >= min_score_chemical
        ) |>
        tidytable::filter(score_final |> as.numeric() >= min_score_final) |>
        keep_best_candidates() |>
        tidytable::inner_join(
          x |>
            load_annotations(show_example = show_example) |>
            tidytable::distinct(
              feature_id,
              inchikey_2D = candidate_structure_inchikey_connectivity_layer
            ) |>
            tidytable::mutate(tidytable::across(
              .cols = tidytable::everything(),
              .fns = function(x) {
                x |> gsub(pattern = "\\|.*", replacement = "")
              }
            ))
        ) |>
        make_confident(score = min_score_final)

      if (predicted_classes) {
        table <- table |>
          tidytable::mutate(
            best_candidate_1 = consensus_1,
            best_candidate_2 = consensus_2,
            best_candidate_3 = consensus_3
          )
      }

      table |>
        tidytable::arrange(tidytable::desc(score_final)) |>
        # ## COMMENT: Dirty way to reduce the number of features
        # ## while keeping some unannotated ones
        # tidytable::mutate(
        #   feature_id = tidytable::coalesce(
        #     inchikey_2D,
        #     rt |>
        #       as.numeric() |>
        #       round(digits = 1L) |>
        #       as.character()
        #   )
        # ) |>
        tidytable::distinct(feature_id, .keep_all = TRUE) |>
        tidytable::mutate(
          intensity = 1,
          comparison_score = 1,
          id = feature_id,
          peak_id = feature_id,
          sample = paste(min_score_initial, min_score_final, name, sep = "_"),
          species = paste(min_score_initial, min_score_final, name, sep = "_"),
          organism = paste(min_score_initial, min_score_final, name, sep = "_")
        )
    })

  return(tables_prepared)
}
