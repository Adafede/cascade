#' Prepare TIMA annotations
#'
#' @export
#'
#' @include keep_best_candidates.R
#' @include load_annotations.R
#' @include make_confident.R
#'
#' @param annotations annotations
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
  tables_prepared <- annotations |>
    purrr::imap(.f = function(x, name) {
      x |>
        load_annotations(show_example = show_example) |>
        tidytable::mutate(tidytable::across(
          .cols = tidyselect::everything(),
          .fns = function(x) {
            x |> gsub(pattern = "\\|.*", replacement = "")
          }
        )) |>
        tidytable::filter(score_initial |> as.numeric() >= min_score_initial) |>
        tidytable::filter(
          score_biological |> as.numeric() >= min_score_biological
        ) |>
        tidytable::filter(
          score_chemical |> as.numeric() >= min_score_chemical
        ) |>
        tidytable::filter(score_final |> as.numeric() >= min_score_final) |>
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
        tidytable::filter(candidate_library %in% libraries) |>
        keep_best_candidates() |>
        tidytable::inner_join(
          x |>
            load_annotations(show_example = show_example) |>
            tidytable::distinct(
              feature_id,
              inchikey_2D = candidate_structure_inchikey_connectivity_layer
            ) |>
            tidytable::mutate(tidytable::across(
              .cols = tidyselect::everything(),
              .fns = function(x) {
                x |> gsub(pattern = "\\|.*", replacement = "")
              }
            ))
        ) |>
        make_confident(score = min_score_final) |>
        tidytable::arrange(tidytable::desc(score_final)) |>
        tidytable::distinct(inchikey_2D, .keep_all = TRUE) |>
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
