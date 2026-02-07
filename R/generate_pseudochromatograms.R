#' Generate pseudochromatograms
#'
#' @export
#'
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
#' @param headers Headers
#' @param detector Detector
#' @param show_example Show example? Default to FALSE
#' @param min_confidence Min confidence
#' @param min_similarity_prefilter Min similarity pre filter
#' @param min_similarity_filter Min similarity filter
#' @param mode Mode
#' @param organism Organism
#' @param fourier_components Fourier components
#' @param frequency Frequency
#' @param resample Resample
#' @param shift Shift
#' @param time_min Time min
#' @param time_max Time max
#' @param intensity_offset Offset to add to intensity values to handle negative
#'   intensities. Default is 100.
#' @param intensity_floor Small value subtracted from minimum intensity.
#'   Default is 0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#' @param baseline_method Method for baseline correction. Default is
#'   "peakDetection".
#'
#' @return A list of plots
#'
#' @examples
#' \dontrun{
#' generate_pseudochromatograms(show_example = TRUE)
#' }
generate_pseudochromatograms <- function(
  annotations = NULL,
  features_informed = NULL,
  features_not_informed = NULL,
  file = NULL,
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  detector = "cad",
  show_example = FALSE,
  min_confidence = 0.4,
  min_similarity_prefilter = 0.6,
  min_similarity_filter = 0.8,
  mode = "pos",
  organism = "Swertia chirayita",
  fourier_components = 0.01,
  frequency = 1,
  resample = 1,
  shift = 0.05,
  time_min = 0.5,
  time_max = 32.5,
  intensity_offset = 100,
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8,
  baseline_method = "peakDetection"
) {
  message("loading annotations")
  annotation_table <- annotations |>
    load_annotations(show_example = show_example)

  message("loading chromatograms")
  chromatograms_all <- file |>
    load_chromatograms(show_example = show_example, headers = headers)

  message("loading name")
  name <- file |>
    load_name(show_example = show_example)

  message("preprocessing chromatograms")
  switch <- switch(
    detector,
    "bpi" = headers["bpi"],
    "cad" = headers["cad"],
    "pda" = headers["pda"]
  )
  list <- chromatograms_all[switch |> names()]
  chromatograms_list <- preprocess_chromatograms(
    detector = detector,
    name = name,
    list = list,
    # signal_name = signal_name,
    shift = shift,
    fourier_components = fourier_components,
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample,
    intensity_offset = intensity_offset,
    intensity_floor = intensity_floor,
    k2 = k2,
    k4 = k4,
    sigma = sigma,
    smoothing_width = smoothing_width,
    baseline_method = baseline_method
  )

  message("keeping only best candidates above desired threshold")
  candidates_confident <- annotation_table |>
    keep_best_candidates() |>
    tidytable::mutate(species = organism) |>
    tidytable::mutate(feature_id = as.numeric(feature_id)) |>
    make_confident(score = min_confidence)

  message("loading compared peaks")
  compared_peaks_list <- prepare_comparison(
    features_informed = features_informed,
    features_not_informed = features_not_informed,
    candidates_confident = candidates_confident,
    min_similarity_prefilter = min_similarity_prefilter,
    min_similarity_filter = min_similarity_filter,
    mode = mode,
    show_example = show_example
  )

  plots_1 <- plot_results_1(
    list = compared_peaks_list,
    chromatogram = chromatograms_list$chromatograms_baselined_long,
    mode = mode,
    time_min = time_min,
    time_max = time_max
  )
  # plots_2 <- plot_results_2(list = compared_peaks_list)

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
      tidytable::bind_rows(
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
    treemaps_progress_no_title(
      xs = names(hierarchies)[
        !grepl(pattern = "_grouped", x = names(hierarchies))
      ],
      hierarchies = hierarchies
    )

  return(list(
    plots_1 = plots_1, # plots_2 = plots_2,
    treemaps = treemaps
  ))
}
