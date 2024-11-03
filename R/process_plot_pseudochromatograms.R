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
                                             resample = 1,
                                             shift = 0.05,
                                             time_min = 0.5,
                                             time_max = 32.5) {
  ANNOTATIONS <- "~/.tima/data/processed/241026_103144_extract/extract_results.tsv"
  BPI <- FALSE
  CAD <- TRUE
  PDA <- FALSE
  CAD_SHIFT <- shift
  PDA_SHIFT <- shift
  TIME_MIN <- time_min
  TIME_MAX <- time_max
  CONFIDENCE_SCORE_MIN <- min_confidence
  PEAK_SIMILARITY_PREFILTER <- min_similarity_prefilter
  PEAK_SIMILARITY <- min_similarity_filter
  IMPORT_FILE_CAD <- "data/interim/peaks/210619_AR_06_V_03_2_01_featuresInformed_cad.tsv"
  IMPORT_FILE_CAD_2 <- "data/interim/peaks/210619_AR_06_V_03_2_01_featuresNotInformed_cad.tsv"
  EXPORT_DIR <- export_dir
  FILE_POSITIVE <- "data/source/mzml/210619_AR_06_V_03_2_01.mzML"
  MODE <- "pos"
  name <- FILE_POSITIVE |>
    gsub(pattern = ".*/", replacement = "") |>
    gsub(pattern = "[0-9]{8}_AR_[0-9]{2}_", replacement = "") |>
    gsub(
      pattern = ".mzML",
      replacement = "",
      fixed = TRUE
    )

  message("listing files")
  message("loading raw files (can take long if loading multiple files)")
  dda_data <- MSnbase::readMSData(
    files = FILE_POSITIVE,
    mode = "onDisk",
    msLevel. = 1
  )

  message("opening raw files objects and extracting chromatograms")
  chromatograms_all <- lapply(FILE_POSITIVE, mzR::openMSfile) |>
    lapply(mzR::chromatograms) |>
    purrr::flatten()


  chromatograms_list_bpi <- preprocess_chromatograms(
    detector = "bpi",
    list = chromatograms_all[c(TRUE, FALSE, FALSE)],
    name = name,
    signal_name = "BasePeak_0",
    shift = 0
  )
  chromatograms_list_cad <- preprocess_chromatograms(
    detector = "cad",
    list = chromatograms_all[c(FALSE, FALSE, TRUE)],
    name = name,
    signal_name = "UV.1_CAD_1_0",
    shift = CAD_SHIFT
  )
  chromatograms_list_pda <-
    preprocess_chromatograms(
      detector = "pda",
      list = chromatograms_all[c(FALSE, TRUE, FALSE)],
      name = name,
      signal_name = "PDA.1_TotalAbsorbance_0",
      shift = PDA_SHIFT
    )


  message("loading annotations")
  annotations <-
    ANNOTATIONS |>
    lapply(
      FUN = function(x) {
        readr::read_delim(file = x) |>
          dplyr::mutate(dplyr::across(
            dplyr::contains("candidate_structure_tax"),
            .fns = function(x) {
              gsub(
                pattern = "\\$",
                replacement = "or",
                x = x
              )
            }
          )) |>
          dplyr::mutate(dplyr::across(
            dplyr::contains("candidate_structure_tax"),
            .fns = function(x) {
              gsub(
                pattern = "u00A7",
                replacement = "$",
                x = x
              )
            }
          )) |>
          dplyr::mutate(mode = MODE)
      }
    ) |>
    dplyr::bind_rows()

  message("keeping best annotations only")
  best_candidates <- annotations |>
    keep_best_candidates()

  message("adding metadata dirtily for now")
  candidates_metadata <- best_candidates |>
    dplyr::mutate(species = "Swertia chirayita") |>
    dplyr::mutate(feature_id = as.numeric(feature_id))

  message("keeping only candidates above desired threshold")
  candidates_confident <- candidates_metadata |>
    make_confident(score = CONFIDENCE_SCORE_MIN)

  myDirtyListF <- function(list, mode = "pos") {
    purrr::map(.x = list, .f = ~ dplyr::filter(.x, grepl(
      pattern = paste0("_", !!as.character(mode)),
      x = sample,
      ignore.case = TRUE
    )))
  }

  compared_peaks_list <- prepare_comparison(
    features_informed = IMPORT_FILE_CAD,
    features_not_informed = IMPORT_FILE_CAD_2,
    min_similarity_prefilter = min_similarity_prefilter,
    min_similarity_filter = min_similarity_filter,
    mode = mode
  )

  plots_1 <- plot_results_1(list = compared_peaks_list)
  plots_2 <- plot_results_2(list = compared_peaks_list)

  fig_taxo <- plots_1$histograms_taxo_maj
  fig_taxo

  fig_minmaj <-
    ggpubr::ggarrange(
      plots_1$histograms_unique_conf_maj,
      ncol = 2,
      nrow = 1,
      plots_1$histograms_unique_conf_min,
      align = "hv",
      common.legend = TRUE,
      legend = "bottom",
      labels = "AUTO"
    )
  fig_minmaj

  hierarchies <- list()
  hierarchies$peaks_maj <-
    prepare_hierarchy(
      dataframe = compared_peaks_list$peaks_maj_precor_taxo_cor |>
        # dplyr::filter(id == "bitter/191113_AR_10043_Pos") |>
        # dplyr::filter(id == "bitter/TODO") |>
        # dplyr::filter(id == "UHR/191109_AR_10043_enriched_UHR_Pos") |>
        dplyr::mutate(id = "peaks_maj") |>
        dplyr::mutate(sample = id, species = id) |>
        dplyr::distinct() |>
        make_other() |>
        no_other(),
      type = "analysis",
      detector = "cad"
    )
  hierarchies$peaks_min <-
    prepare_hierarchy(
      dataframe = compared_peaks_list$peaks_min_precor_taxo_cor |>
        dplyr::filter(mode == "pos") |>
        dplyr::mutate(id = "peaks_min") |>
        dplyr::mutate(sample = id, species = id) |>
        dplyr::distinct() |>
        make_other() |>
        no_other() |>
        dplyr::mutate(intensity = 1),
      # because MS intensity not OK
      type = "analysis",
      detector = "ms"
    )
  hierarchies$special <- prepare_hierarchy(
    dataframe = compared_peaks_list$peaks_maj_precor_taxo_cor |>
      # dplyr::filter(id == "bitter/191113_AR_10043_Pos") |>
      # dplyr::filter(id == "bitter/TODO") |>
      # dplyr::filter(id == "UHR/191109_AR_10043_enriched_UHR_Pos") |>
      dplyr::filter(mode == "pos") |>
      dplyr::mutate(id = "peaks_maj") |>
      dplyr::mutate(sample = id, species = id) |>
      dplyr::select(-taxo, -taxo_2, -sum, -sum_2, -keep) |>
      rbind(
        compared_peaks_list$peaks_min_precor_taxo_cor |>
          # dplyr::filter(id == "UHR/191109_AR_10043_enriched_UHR_Pos") |>
          dplyr::filter(mode == "pos") |>
          dplyr::mutate(id = "peaks_min") |>
          dplyr::mutate(sample = id, species = id)
      ) |>
      dplyr::distinct(),
    type = "analysis",
    detector = "cad"
  ) |>
    dplyr::arrange(sample)
  treemaps <-
    treemaps_progress_no_title(xs = names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))], hierarchies = hierarchies)
  # treemaps$special
}
