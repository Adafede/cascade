#' Process plot pseudochromatograms
#'
#' @export
#'
#' @include check_export_dir.R
#' @include keep_best_candidates.R
#' @include make_confident.R
#' @include plot_results.R
#' @include prepare_comparison.R
#' @include prepare_hierarchy.R
#' @include preprocess_chromatograms.R
#' @include treemaps_progress.R
#'
#' @param annotations Annotations
#' @param features_informed Features informed
#' @param features_not_informed Features not informed
#' @param file File
#' @param detector Detector
#' @param export_dir Export dir
#' @param show_example Show example? Default to FALSE
#' @param min_confidence Min confidence
#' @param min_similarity_prefilter Min similarity pre filter
#' @param min_similarity_filter Min similarity filter
#' @param mode Mode
#' @param organism Organism
#' @param resample Resample
#' @param shift Shift
#' @param time_min Time min
#' @param time_max Time max
#'
#' @return A plot with histograms on chromatograms
#'
#' @examples
#' \dontrun{
#' process_plot_pseudochromatograms(show_example = TRUE)
#' }
process_plot_pseudochromatograms <- function(annotations = NULL,
                                             features_informed = NULL,
                                             features_not_informed = NULL,
                                             file = NULL,
                                             detector = "cad",
                                             export_dir = "data/figures",
                                             show_example = FALSE,
                                             min_confidence = 0.4,
                                             min_similarity_prefilter = 0.6,
                                             min_similarity_filter = 0.8,
                                             mode = "pos",
                                             organism = "Swertia chirayita",
                                             resample = 1,
                                             shift = 0.05,
                                             time_min = 0.5,
                                             time_max = 32.5) {
  EXPORT_DIR <- export_dir

  message("loading compared peaks")
  compared_peaks_list <- prepare_comparison(
    features_informed = features_informed,
    features_not_informed = features_not_informed,
    min_similarity_prefilter = min_similarity_prefilter,
    min_similarity_filter = min_similarity_filter,
    mode = mode,
    show_example = show_example
  )

  message("loading annotations")
  annotation_table <- annotations |>
    load_annotations(show_example = show_example)

  message("keeping best annotations only")
  best_candidates <- annotation_table |>
    keep_best_candidates()

  message("adding metadata")
  candidates_metadata <- best_candidates |>
    dplyr::mutate(species = organism) |>
    dplyr::mutate(feature_id = as.numeric(feature_id))

  message("keeping only candidates above desired threshold")
  candidates_confident <- candidates_metadata |>
    make_confident(score = min_confidence)

  myDirtyListF <- function(list, mode = "pos") {
    purrr::map(.x = list, .f = ~ dplyr::filter(.x, grepl(
      pattern = paste0("_", !!as.character(mode)),
      x = sample,
      ignore.case = TRUE
    )))
  }

  plots_1 <- plot_results_1(list = compared_peaks_list)
  plots_2 <- plot_results_2(list = compared_peaks_list)

  fig_taxo <- plots_1$histograms_taxo_maj

  hierarchies <- list()
  hierarchies$peaks_maj <-
    prepare_hierarchy(
      dataframe = compared_peaks_list$peaks_maj_precor_taxo_cor |>
        tidytable::mutate(id = "peaks_maj") |>
        tidytable::mutate(sample = id, species = id) |>
        tidytable::distinct() |>
        make_other() |>
        no_other(),
      type = "analysis",
      detector = detector
    )
  hierarchies$peaks_min <-
    prepare_hierarchy(
      dataframe = compared_peaks_list$peaks_min_precor_taxo_cor |>
        tidytable::filter(mode == "pos") |>
        tidytable::mutate(id = "peaks_min") |>
        tidytable::mutate(sample = id, species = id) |>
        tidytable::distinct() |>
        make_other() |>
        no_other() |>
        tidytable::mutate(intensity = 1),
      # because MS intensity not OK
      type = "analysis",
      detector = "ms"
    )
  hierarchies$special <- prepare_hierarchy(
    dataframe = compared_peaks_list$peaks_maj_precor_taxo_cor |>
      tidytable::filter(mode == mode) |>
      tidytable::mutate(id = "peaks_maj") |>
      tidytable::mutate(sample = id, species = id) |>
      tidytable::select(-taxo, -taxo_2, -sum, -sum_2, -keep) |>
      rbind(
        compared_peaks_list$peaks_min_precor_taxo_cor |>
          tidytable::filter(mode == mode) |>
          tidytable::mutate(id = "peaks_min") |>
          tidytable::mutate(sample = id, species = id)
      ) |>
      tidytable::distinct(),
    type = "analysis",
    detector = detector
  ) |>
    tidytable::arrange(sample)

  treemaps <-
    treemaps_progress_no_title(xs = names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))], hierarchies = hierarchies)
}
