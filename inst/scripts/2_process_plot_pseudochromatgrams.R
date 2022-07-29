start <- Sys.time()

library(package = data.table, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = ggalluvial, quietly = TRUE)
library(package = ggfittext, quietly = TRUE)
library(package = ggplot2, quietly = TRUE)
library(package = ggpubr, quietly = TRUE)
library(package = parallel, quietly = TRUE)
library(package = patchwork, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = readr, quietly = TRUE)

source(file = "R/add_peak_metadata.R")
source(file = "R/check_export_dir.R")
source(file = "R/colors.R")
source(file = "R/get_gnps.R")
source(file = "R/get_params.R")
source(file = "R/keep_best_candidates.R")
source(file = "R/log_debug.R")
source(file = "R/make_confident.R")
source(file = "R/make_other.R")
source(file = "R/no_other.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/plot_histograms.R")
source(file = "R/plot_peaks_statistics.R")
source(file = "R/plot_results.R")
source(file = "R/prehistograms_progress.R")
source(file = "R/prepare_comparison.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_hierarchy_preparation.R")
source(file = "R/prepare_plot.R")
source(file = "R/treemaps_progress.R")
source(file = "R/y_as_na.R")

step <- "processing"
paths <- parse_yaml_paths()
params <- ""
params <- get_params(step = step)

log_debug(
  "This program performs",
  "TODO"
)
log_debug("Authors: \n", "AR")
log_debug("Contributors: \n", "...")

#' Dirty generic paths and parameters
source(file = "R/dirty_paths_params.R")

#' Specific paths
TIME_MIN <- 0.5
TIME_MAX <- 127


###
#' DIRTY FROM 1_process_compare_peaks.R
source(file = "R/baseline_chromatogram.R")
source(file = "R/change_intensity_name.R")
source(file = "R/check_export_dir.R")
source(file = "R/compare_peaks.R")
source(file = "R/extract_ms.R")
source(file = "R/extract_ms_peak.R")
source(file = "R/filter_ms.R")
source(file = "R/get_gnps.R")
source(file = "R/get_params.R")
source(file = "R/improve_signal.R")
source(file = "R/improve_signals_progress.R")
source(file = "R/join_peaks.R")
source(file = "R/log_debug.R")
source(file = "R/make_confident.R")
source(file = "R/normalize_chromato.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/peaks_progress.R")
source(file = "R/plot_chromatogram.R")
source(file = "R/prepare_features.R")
source(file = "R/prepare_mz.R")
source(file = "R/prepare_peaks.R")
source(file = "R/prepare_rt.R")
source(file = "R/preprocess_chromatograms.R")
source(file = "R/preprocess_peaks.R")
source(file = "R/process_peaks.R")
source(file = "R/transform_ms.R")
source(file = "R/y_as_na.R")
source(file = "R/dirty_paths_params.R")
log_debug(x = "listing files")
files <- list.files(
  path = TOYSET,
  pattern = paste0(params$filename$mzml, ".mzML"),
  full.names = TRUE,
  recursive = TRUE
)
names <- list.files(
  path = TOYSET,
  pattern = paste0(params$filename$mzml, ".mzML"),
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML",
    replacement = "",
    fixed = TRUE
  )
log_debug(x = "loading raw files (can take long if loading multiple files)")
dda_data <- MSnbase::readMSData(
  files = files,
  mode = "onDisk",
  msLevel. = 1
)
log_debug(x = "opening raw files objects and extracting chromatograms")
chromatograms_all <- lapply(files, mzR::openMSfile) |>
  lapply(mzR::chromatograms) |>
  purrr::flatten()
chromatograms_list_bpi <- preprocess_chromatograms(
  detector = "bpi",
  list = chromatograms_all[c(TRUE, FALSE, FALSE)],
  signal_name = "BasePeak_0",
  shift = 0
)
chromatograms_list_cad <- preprocess_chromatograms()
chromatograms_list_pda <-
  preprocess_chromatograms(
    detector = "pda",
    list = chromatograms_all[c(FALSE, TRUE, FALSE)],
    signal_name = "PDA.1_TotalAbsorbance_0",
    shift = PDA_SHIFT
  )
###
log_debug(x = "loading annotations")
#' TODO add ionization mode?
annotations <-
  ANNOTATIONS |>
  lapply(
    FUN = function(x) {
      readr::read_delim(file = x) |>
        dplyr::mutate(best_candidate = gsub(
          pattern = "§",
          replacement = "$",
          x = best_candidate
        )) |>
        dplyr::mutate(mode = ifelse(
          test = grepl(
            pattern = "_pos.",
            x = x,
            ignore.case = TRUE
          ),
          yes = "pos",
          no = "neg"
        ))
    }
  ) |>
  dplyr::bind_rows()

log_debug(x = "keeping best annotations only")
best_candidates <- annotations |>
  keep_best_candidates()

log_debug(x = "adding metadata dirtily for now")
candidates_metadata <- best_candidates |>
  dplyr::mutate(species = "Swertia chirayita") |>
  dplyr::mutate(feature_id = as.numeric(feature_id))

log_debug(x = "keeping only candidates above desired threshold")
candidates_confident <- candidates_metadata |>
  make_confident(score = CONFIDENCE_SCORE_MIN)

myDirtyListF <- function(list, mode = "pos") {
  purrr::map(.x = list, .f = ~ dplyr::filter(.x, grepl(
    pattern = paste0("_", !!as.character(mode)),
    x = sample,
    ignore.case = TRUE
  )))
}

if (params$signal$detector$bpi == TRUE) {
  detector <- "bpi"
}
if (params$signal$detector$bpi == TRUE) {
  compared_peaks_list_bpi <- prepare_comparison(detector = "bpi")
  compared_peaks_list_bpi_pos <- compared_peaks_list_bpi |>
    myDirtyListF()
  compared_peaks_list_bpi_neg <- compared_peaks_list_bpi |>
    myDirtyListF(mode = "neg")
  plots_1_bpi_pos <- plot_results_1(detector = "bpi")
  plots_1_bpi_neg <- plot_results_1(detector = "bpi", mode = "neg")
  plots_2_bpi <- plot_results_2(detector = "bpi")
}
if (params$signal$detector$cad == TRUE) {
  detector <- "cad"
}
if (params$signal$detector$cad == TRUE) {
  compared_peaks_list_cad <- prepare_comparison()
  compared_peaks_list_cad_pos <- compared_peaks_list_cad |>
    myDirtyListF()
  compared_peaks_list_cad_neg <- compared_peaks_list_cad |>
    myDirtyListF(mode = "neg")
  plots_1_cad_pos <- plot_results_1()
  plots_1_cad_neg <- plot_results_1(mode = "neg")
  plots_2_cad <- plot_results_2()
}

if (params$signal$detector$pda == TRUE) {
  detector <- "pda"
}

if (params$signal$detector$pda == TRUE) {
  compared_peaks_list_pda <- prepare_comparison(detector = "pda")
  compared_peaks_list_pda_pos <- compared_peaks_list_pda |>
    myDirtyListF()
  compared_peaks_list_pda_neg <- compared_peaks_list_pda |>
    myDirtyListF(mode = "neg")
  plots_1_pda_pos <- plot_results_1(detector = "pda")
  plots_1_pda_neg <- plot_results_1(detector = "pda", mode = "neg")
  plots_2_pda <- plot_results_2(detector = "pda")
}

fig_taxo <-
  ggpubr::ggarrange(
    plots_1_cad_pos$histograms_taxo_maj,
    ncol = 1,
    nrow = 2,
    plots_1_cad_neg$histograms_taxo_maj,
    align = "hv",
    common.legend = TRUE,
    legend = "bottom",
    labels = "AUTO"
  )
fig_taxo

fig_minmaj <-
  ggpubr::ggarrange(
    plots_1_cad_pos$histograms_unique_conf_maj,
    ncol = 2,
    nrow = 1,
    plots_1_cad_pos$histograms_unique_conf_min,
    align = "hv",
    # common.legend = TRUE,
    legend = "bottom"
  )
fig_minmaj

# test_taxo_3 <-
#   ggpubr::ggarrange(
#     ncol = 2,
#     nrow = 2,
#     plots_1_cad_pos$histograms_unique_conf_maj,
#     plots_1_cad_neg$histograms_unique_conf_maj,
#     plots_1_cad_pos$histograms_unique_conf_min,
#     plots_1_cad_neg$histograms_unique_conf_min,
#     common.legend = TRUE,
#     legend = "bottom",
#     labels = "AUTO"
#   )
# test_taxo_3

hierarchies <- list()
hierarchies$peaks_maj <-
  prepare_hierarchy(
    dataframe = compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
      dplyr::filter(id == "UHR/191109_AR_10043_enriched_UHR_Pos") |>
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
    dataframe = compared_peaks_list_cad$peaks_min_precor_taxo_cor |>
      dplyr::filter(mode == "pos") |>
      dplyr::mutate(id = "peaks_min") |>
      dplyr::mutate(sample = id, species = id) |>
      dplyr::distinct() |>
      make_other() |>
      no_other(),
    type = "analysis",
    detector = "ms"
  )
hierarchies$special <- prepare_hierarchy(
  dataframe = compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
    dplyr::filter(id == "UHR/191109_AR_10043_enriched_UHR_Pos") |>
    dplyr::mutate(id = "peaks_maj") |>
    dplyr::mutate(sample = id, species = id) |>
    dplyr::select(-taxo, -taxo_2, -sum, -sum_2, -keep) |>
    rbind(
      compared_peaks_list_cad$peaks_min_precor_taxo_cor |>
        dplyr::filter(mode == "pos") |>
        dplyr::mutate(id = "peaks_min") |>
        dplyr::mutate(sample = id, species = id)
    ) |>
    dplyr::distinct(),
  type = "analysis",
  detector = "cad"
)
treemaps <-
  treemaps_progress_noTitle(xs = names(hierarchies)[!grepl(
    pattern = "_grouped",
    x = names(hierarchies)
  )])
treemaps$special

df_meta_bpi_pos <- compared_peaks_list_bpi$peaks_all |>
  dplyr::filter(mode == "pos") |>
  dplyr::full_join(
    best_candidates |>
      dplyr::filter(mode == "pos") |>
      dplyr::mutate(feature_id = as.numeric(feature_id))
  ) |>
  dplyr::left_join(
    compared_peaks_list_bpi$peaks_maj_precor_taxo_cor |>
      dplyr::distinct(sample, peak_id, mode, feature_id, keep)
  ) |>
  add_peak_metadata()
df_meta_bpi_neg <- compared_peaks_list_bpi$peaks_all |>
  dplyr::filter(mode == "neg") |>
  dplyr::full_join(
    best_candidates |>
      dplyr::filter(mode == "neg") |>
      dplyr::mutate(feature_id = as.numeric(feature_id))
  ) |>
  dplyr::left_join(
    compared_peaks_list_bpi$peaks_maj_precor_taxo_cor |>
      dplyr::distinct(sample, peak_id, mode, feature_id, keep)
  ) |>
  add_peak_metadata()

df_meta_cad_pos <- compared_peaks_list_cad$peaks_all |>
  dplyr::filter(mode == "pos") |>
  dplyr::full_join(
    best_candidates |>
      dplyr::filter(mode == "pos") |>
      dplyr::mutate(feature_id = as.numeric(feature_id))
  ) |>
  dplyr::left_join(
    compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
      dplyr::distinct(sample, peak_id, mode, feature_id, keep)
  ) |>
  add_peak_metadata()
df_meta_cad_neg <- compared_peaks_list_cad$peaks_all |>
  dplyr::filter(mode == "neg") |>
  dplyr::full_join(
    best_candidates |>
      dplyr::filter(mode == "neg") |>
      dplyr::mutate(feature_id = as.numeric(feature_id))
  ) |>
  dplyr::left_join(
    compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
      dplyr::distinct(sample, peak_id, mode, feature_id, keep)
  ) |>
  add_peak_metadata()

df_meta_pda_pos <- compared_peaks_list_pda$peaks_all |>
  dplyr::filter(mode == "pos") |>
  dplyr::full_join(
    best_candidates |>
      dplyr::filter(mode == "pos") |>
      dplyr::mutate(feature_id = as.numeric(feature_id))
  ) |>
  dplyr::left_join(
    compared_peaks_list_pda$peaks_maj_precor_taxo_cor |>
      dplyr::distinct(sample, peak_id, mode, feature_id, keep)
  ) |>
  add_peak_metadata()
df_meta_pda_neg <- compared_peaks_list_pda$peaks_all |>
  dplyr::filter(mode == "neg") |>
  dplyr::full_join(
    best_candidates |>
      dplyr::filter(mode == "neg") |>
      dplyr::mutate(feature_id = as.numeric(feature_id))
  ) |>
  dplyr::left_join(
    compared_peaks_list_pda$peaks_maj_precor_taxo_cor |>
      dplyr::distinct(sample, peak_id, mode, feature_id, keep)
  ) |>
  add_peak_metadata()

plots_bpi_pos <- df_meta_bpi_pos |>
  plot_peaks_statistics()
plots_bpi_neg <- df_meta_bpi_neg |>
  plot_peaks_statistics()

plots_cad_pos <- df_meta_cad_pos |>
  plot_peaks_statistics()
plots_cad_neg <- df_meta_cad_neg |>
  plot_peaks_statistics()

plots_pda_pos <- df_meta_pda_pos |>
  plot_peaks_statistics()
plots_pda_neg <- df_meta_pda_neg |>
  plot_peaks_statistics()

#' export
# ggplot2::ggsave(
#   filename = "~/git/cascade/data/paper/cascade-5.pdf",
#   plot = fig_taxo
# )
# ggplot2::ggsave(
#   filename = "~/git/cascade/data/paper/cascade-6.pdf",
#   plot = ggpubr::ggarrange(
#     plots_cad_pos[[1]],
#     plots_cad_pos[[3]],
#     labels = "AUTO",
#     nrow = 2,
#     common.legend = FALSE,
#     legend = "bottom"
#   ),
#   width = 8,
#   height = 9,
#   limitsize = FALSE
# )
# ggplot2::ggsave(
#   filename = "~/git/cascade/data/paper/cascade-7-a.pdf",
#   plot = fig_minmaj,
#   width = 32,
#   height = 18,
#   limitsize = FALSE
# )
# plotly::save_image(
#   p = treemaps$special,
#   file = "data/paper/cascade-7-b.pdf",
#   width = 1600,
#   height = 900
# )

#' WIP
readr::write_csv(
  x = compared_peaks_list_cad_pos[["peaks_maj_precor_taxo_cor"]] |>
    dplyr::select(
      peak_id,
      peak_area,
      comparison_score,
      feature_id
    ),
  file = "inst/extdata/interim/peaks/191109_AR_10043_enriched_UHR_Pos_featuresInformed_filtered_cad.csv"
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
