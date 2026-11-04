#' Generate annotated pseudochromatograms and chemical class treemaps
#'
#' Combines preprocessed detector chromatograms, peak-correlated MS features,
#' and TIMA structural annotations to produce a suite of annotated
#' pseudochromatograms (colored by chemical pathway/superclass and closest
#' taxonomic occurrence) and interactive chemical hierarchy treemaps.
#'
#' @export
#'
#' @param annotations Path to TIMA annotations `.tsv` file (optional when
#'   `show_example = TRUE`).
#' @param features_informed Path to peak-informed features `.tsv` file
#'   (optional when `show_example = TRUE`).
#' @param features_not_informed Path to non-informed features `.tsv` file
#'   (optional when `show_example = TRUE`).
#' @param file Path to `.mzML` file (optional when `show_example = TRUE`).
#' @param headers Named character vector mapping detector types to `.mzML`
#'   chromatogram header IDs.
#' @param detector Detector type: `"cad"`, `"bpi"`, or `"pda"`. Default is
#'   `"cad"`.
#' @param show_example Logical; use bundled example data. Default is `FALSE`.
#' @param min_confidence Minimum final annotation confidence score. Default is
#'   `0.4`.
#' @param min_similarity_prefilter Minimum correlation score for pre-filtering.
#'   Default is `0.6`.
#' @param min_similarity_filter Minimum correlation score for final filtering.
#'   Default is `0.8`.
#' @param mode Ionization mode: `"pos"` or `"neg"`. Default is `"pos"`.
#' @param organism Investigated organism scientific name for taxonomic
#'   prioritization. Default is `"Swertia chirayita"`.
#' @param fourier_components Fraction of Fourier components to retain. Default
#'   is `0.01`.
#' @param frequency Acquisition frequency in Hz. Default is `1`.
#' @param resample Resampling factor. Default is `1`.
#' @param shift Retention time shift in minutes. Default is `0.05`.
#' @param time_min Minimum retention time in minutes. Default is `0.5`.
#' @param time_max Maximum retention time in minutes. Default is `32.5`.
#' @param intensity_floor Small positive value for intensity floor. Default is
#'   `0.001`.
#' @param k2 Second-derivative weight for signal sharpening. Default is `250`.
#' @param k4 Fourth-derivative weight for signal sharpening. Default is
#'   `1250000`.
#' @param sigma Gaussian sigma for signal sharpening. Default is `0.05`.
#' @param smoothing_width Moving-average window width for derivative smoothing.
#'   Default is `8`.
#' @param baseline_method Baseline correction method passed to
#'   [baseline::baseline()]. Default is `"peakDetection"`.
#' @param improve_signal Logical; apply Fourier filtering and sharpening.
#'   Default is `TRUE`.
#' @param correct_acn Logical; apply acetonitrile gradient response correction
#'   to CAD peak areas. Default is `FALSE`.
#' @param acn_time Numeric; gradient time in minutes to reach 100% acetonitrile
#'   when `correct_acn = TRUE`. Default is `30.5`.
#' @param dvol Numeric; dwell time offset in minutes when `correct_acn = TRUE`.
#'   Default is `0.5`.
#'
#' @return A named list with three elements:
#'   - `plots_1`: Named list of 8 `ggplot` pseudochromatograms.
#'   - `plots_2`: Named list of 2 interactive `plotly` signal-based sunbursts.
#'   - `treemaps`: Named list of interactive `plotly` treemap widgets.
#'
#' @concept EDAM:operation_3203 Chromatogram visualisation
#' @concept EDAM:operation_2942 Treemap visualisation
#' @concept EDAM:operation_3803 Natural product identification
#' @concept EDAM:topic_3172 Metabolomics
#' @concept EDAM:topic_0154 Small molecules
#'
#' @examples
#' \donttest{
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
  intensity_floor = 0.001,
  k2 = 250,
  k4 = 1250000,
  sigma = 0.05,
  smoothing_width = 8,
  baseline_method = "peakDetection",
  improve_signal = TRUE,
  correct_acn = FALSE,
  acn_time = 30.5,
  dvol = 0.5
) {
  target_mode <- mode

  message("loading annotations")
  annotation_table <- load_annotations(annotations, show_example = show_example)

  message("loading chromatograms")
  chromatograms_all <- load_chromatograms(
    file,
    show_example = show_example,
    headers = headers
  )

  message("loading name")
  name <- load_name(file, show_example = show_example)

  message("preprocessing chromatograms")
  selected_header <- switch(
    detector,
    "bpi" = headers["bpi"],
    "cad" = headers["cad"],
    "pda" = headers["pda"]
  )
  chrom_sub <- chromatograms_all[names(selected_header)]
  chromatograms_list <- preprocess_chromatograms(
    detector = detector,
    name = name,
    list = chrom_sub,
    shift = shift,
    fourier_components = fourier_components,
    time_min = time_min,
    time_max = time_max,
    frequency = frequency,
    resample = resample,
    intensity_floor = intensity_floor,
    k2 = k2,
    k4 = k4,
    sigma = sigma,
    smoothing_width = smoothing_width,
    baseline_method = baseline_method,
    improve_signal = improve_signal
  )

  message("keeping only best candidates above desired threshold")
  candidates_confident <- keep_best_candidates(annotation_table)
  candidates_confident[, `:=`(
    species = organism,
    feature_id = as.numeric(feature_id)
  )]
  candidates_confident <- make_confident(
    candidates_confident,
    score = min_confidence
  )

  message("loading compared peaks")
  compared_peaks_list <- prepare_comparison(
    features_informed = features_informed,
    features_not_informed = features_not_informed,
    candidates_confident = candidates_confident,
    min_similarity_prefilter = min_similarity_prefilter,
    min_similarity_filter = min_similarity_filter,
    mode = target_mode,
    show_example = show_example,
    correct_acn = correct_acn,
    acn_time = acn_time,
    dvol = dvol
  )

  plots_1 <- plot_results_1(
    list = compared_peaks_list,
    chromatogram = chromatograms_list$chromatograms_baselined_long,
    mode = target_mode,
    time_min = time_min,
    time_max = time_max
  )

  plots_2 <- plot_results_2(list = compared_peaks_list)

  hierarchies <- list()

  maj_dt <- data.table::copy(compared_peaks_list$peaks_maj_precor_taxo_cor)
  maj_dt[, `:=`(id = "peaks_maj", sample = "peaks_maj", species = "peaks_maj")]
  hierarchies$peaks_maj <- prepare_hierarchy(
    dataframe = no_other(make_other(unique(maj_dt))),
    type = "analysis",
    detector = detector
  )

  min_dt <- data.table::copy(
    compared_peaks_list$peaks_min_precor_taxo_cor[mode == target_mode]
  )
  min_dt[, `:=`(id = "peaks_min", sample = "peaks_min", species = "peaks_min")]
  min_dt_clean <- no_other(make_other(unique(min_dt)))
  min_dt_clean[, intensity := 1]
  hierarchies$peaks_min <- prepare_hierarchy(
    dataframe = min_dt_clean,
    type = "analysis",
    detector = "ms"
  )

  maj_spec <- data.table::copy(
    compared_peaks_list$peaks_maj_precor_taxo_cor[mode == target_mode]
  )
  maj_spec[, `:=`(
    id = "peaks_maj",
    sample = "peaks_maj",
    species = "peaks_maj"
  )]
  drop_tax_cols <- intersect(
    c("taxo", "taxo_2", "sum", "sum_2", "keep"),
    names(maj_spec)
  )
  if (length(drop_tax_cols) > 0L) {
    maj_spec[, (drop_tax_cols) := NULL]
  }

  min_spec <- data.table::copy(
    compared_peaks_list$peaks_min_precor_taxo_cor[mode == target_mode]
  )
  min_spec[, `:=`(
    id = "peaks_min",
    sample = "peaks_min",
    species = "peaks_min"
  )]

  spec_combined <- unique(data.table::rbindlist(
    list(maj_spec, min_spec),
    use.names = TRUE,
    fill = TRUE
  ))
  h_special <- prepare_hierarchy(
    dataframe = spec_combined,
    type = "analysis",
    detector = detector
  )
  if (nrow(h_special) > 0L) {
    data.table::setorder(h_special, sample, na.last = TRUE)
  }
  hierarchies$special <- h_special

  treemaps <- treemaps_progress_no_title(
    xs = names(hierarchies)[!grepl("_grouped", names(hierarchies))],
    hierarchies = hierarchies
  )

  list(
    plots_1 = plots_1,
    plots_2 = plots_2,
    treemaps = treemaps
  )
}
