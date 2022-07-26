start <- Sys.time()

library(package = data.table, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
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
source(file = "R/prepare_comparison.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_hierarchy_preparation.R")
source(file = "R/prepare_plot.R")
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

if (params$signal$detector$bpi == TRUE) {
  detector <- "bpi"
}
if (params$signal$detector$bpi == TRUE) {
  compared_peaks_list_bpi <- prepare_comparison(detector = "bpi")
  plots_1_bpi <- plot_results_1(detector = "bpi")
  plots_2_bpi <- plot_results_2(detector = "bpi")
}
if (params$signal$detector$cad == TRUE) {
  detector <- "cad"
}
if (params$signal$detector$cad == TRUE) {
  compared_peaks_list_cad <- prepare_comparison()
  plots_1_cad <- plot_results_1()
  plots_2_cad <- plot_results_2()
}
if (params$signal$detector$pda == TRUE) {
  detector <- "pda"
}
if (params$signal$detector$pda == TRUE) {
  compared_peaks_list_pda <- prepare_comparison(detector = "pda")
  plots_1_pda <- plot_results_1(detector = "pda")
  plots_2_pda <- plot_results_2(detector = "pda")
}

#' Work in progress
test <-
  ggpubr::ggarrange(
    plots_1_cad$histograms_taxo_maj,
    ncol = 1,
    nrow = 2,
    plots_1_cad$histograms_taxo_min,
    align = "hv"
  )

#' Work in progress
#' See how to do best also with non-annotated peaks, etc.
df_meta_bpi_pos <-
  compared_peaks_list_bpi$peaks_maj_precor_taxo_cor |>
  dplyr::filter(grepl(
    pattern = "_pos",
    x = sample,
    ignore.case = TRUE
  )) |>
  add_peak_metadata()
df_meta_bpi_neg <-
  compared_peaks_list_bpi$peaks_maj_precor_taxo_cor |>
  dplyr::filter(grepl(
    pattern = "_neg",
    x = sample,
    ignore.case = TRUE
  )) |>
  add_peak_metadata()

df_meta_cad_pos <-
  compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
  dplyr::filter(grepl(
    pattern = "_pos",
    x = sample,
    ignore.case = TRUE
  )) |>
  add_peak_metadata()
df_meta_cad_neg <-
  compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
  dplyr::filter(grepl(
    pattern = "_neg",
    x = sample,
    ignore.case = TRUE
  )) |>
  add_peak_metadata()

df_meta_pda_pos <-
  compared_peaks_list_pda$peaks_maj_precor_taxo_cor |>
  dplyr::filter(grepl(
    pattern = "_pos",
    x = sample,
    ignore.case = TRUE
  )) |>
  add_peak_metadata()
df_meta_pda_neg <-
  compared_peaks_list_pda$peaks_maj_precor_taxo_cor |>
  dplyr::filter(grepl(
    pattern = "_neg",
    x = sample,
    ignore.case = TRUE
  )) |>
  add_peak_metadata()

df_meta_bpi_pos |>
  plot_peaks_statistics()
df_meta_bpi_neg |>
  plot_peaks_statistics()

df_meta_cad_pos |>
  plot_peaks_statistics()
df_meta_cad_neg |>
  plot_peaks_statistics()

df_meta_pda_pos |>
  plot_peaks_statistics()
df_meta_pda_neg |>
  plot_peaks_statistics()

end <- Sys.time()

log_debug("Script finished in", format(end - start))
