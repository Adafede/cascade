#' Prepare TIMA structural and predicted annotation tables
#'
#' Filters and standardizes Taxonomically Informed Metabolite Annotation (TIMA)
#' output tables across initial spectral similarity, matched peak counts,
#' biological, chemical, and final confidence scores, retaining the top
#' candidate per feature for downstream visualization with [plot_tima()].
#'
#' @export
#'
#' @param annotations Named character vector of annotation file paths (optional
#'   when `show_example = TRUE`).
#' @param predicted_classes Logical; use Sirius/CANOPUS predicted NPClassifier
#'   classes instead of candidate structure classes. Default is `FALSE`.
#' @param min_score_initial Minimum initial spectral similarity score. Default
#'   is `0.0`.
#' @param min_score_biological Minimum biological (taxonomic) score. Default is
#'   `0.0`.
#' @param min_score_chemical Minimum chemical consistency score. Default is
#'   `0.0`.
#' @param min_score_final Minimum final combined score. Default is `0.0`.
#' @param min_matched_peaks_absolute Minimum absolute number of matched MS2
#'   peaks. Default is `0L`.
#' @param min_matched_peaks_percentage Minimum fraction of spectrum peaks
#'   matched. Default is `0.0`.
#' @param min_peaks Minimum total number of peaks in the query MS2 spectrum.
#'   Default is `3L`.
#' @param libraries Character vector of spectral libraries to retain. Default
#'   is `c("gnps", "massbank", "merlin", "ISDB", "ISDB - Wikidata", "TIMA MS1")`.
#' @param show_example Logical; use bundled example data. Default is `FALSE`.
#'
#' @return A named list of prepared `data.table` annotation tables.
#'
#' @concept EDAM:operation_3803 Natural product identification
#' @concept EDAM:operation_3214 Spectral analysis
#' @concept EDAM:operation_3801 Spectral library search
#' @concept EDAM:data_0962 Small molecule report
#' @concept EDAM:topic_3172 Metabolomics
#' @concept EDAM:topic_0154 Small molecules
#'
#' @examples
#' prepare_tima_annotations(show_example = TRUE)
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
  if (is.null(annotations)) {
    annotations <- c("annotations" = tempfile())
  }
  if (predicted_classes) {
    names(annotations) <- paste(names(annotations), "pred", sep = "_")
  }

  nms <- names(annotations)
  if (is.null(nms)) {
    nms <- as.character(seq_along(annotations))
  }

  stats::setNames(
    Map(
      function(x, name) {
        dt <- data.table::copy(data.table::as.data.table(
          load_annotations(x, show_example = show_example)
        ))
        for (col in names(dt)) {
          data.table::set(dt, j = col, value = gsub("\\|.*", "", dt[[col]]))
        }

        if (!predicted_classes) {
          s_init <- suppressWarnings(as.numeric(dt$score_initial))
          peaks_matched <- suppressWarnings(as.numeric(
            dt$candidate_count_similarity_peaks_matched
          ))
          spec_peaks <- suppressWarnings(as.numeric(dt$feature_spectrum_peaks))

          keep_idx <- !is.na(s_init) &
            s_init >= min_score_initial &
            ((!is.na(peaks_matched) &
              !is.na(spec_peaks) &
              peaks_matched >= spec_peaks * min_matched_peaks_percentage) |
              (!is.na(peaks_matched) &
                peaks_matched >= min_matched_peaks_absolute)) &
            !is.na(spec_peaks) &
            spec_peaks >= min_peaks &
            dt$candidate_library %in% libraries
          dt <- dt[keep_idx]
        } else {
          dt[, `:=`(
            candidate_structure_tax_npc_01pat = feature_pred_tax_npc_01pat_val,
            candidate_structure_tax_npc_02sup = feature_pred_tax_npc_02sup_val,
            candidate_structure_tax_npc_03cla = feature_pred_tax_npc_03cla_val
          )]
        }

        s_bio <- suppressWarnings(as.numeric(dt$score_biological))
        s_chem <- suppressWarnings(as.numeric(dt$score_chemical))
        s_fin <- suppressWarnings(as.numeric(dt$score_final))

        dt <- dt[
          !is.na(s_bio) &
            s_bio >= min_score_biological &
            !is.na(s_chem) &
            s_chem >= min_score_chemical &
            !is.na(s_fin) &
            s_fin >= min_score_final
        ]

        table <- make_confident(
          keep_best_candidates(dt),
          score = min_score_final
        )

        if (predicted_classes) {
          table[, `:=`(
            best_candidate_1 = consensus_1,
            best_candidate_2 = consensus_2,
            best_candidate_3 = consensus_3
          )]
        }

        data.table::setorder(table, -score_final, na.last = TRUE)
        table <- unique(table, by = "feature_id")
        sample_label <- paste(
          min_score_initial,
          min_score_final,
          name,
          sep = "_"
        )
        table[, `:=`(
          intensity = 1,
          comparison_score = 1,
          id = feature_id,
          peak_id = feature_id,
          sample = sample_label,
          species = sample_label,
          organism = sample_label
        )]
        table
      },
      annotations,
      nms
    ),
    nms
  )
}
