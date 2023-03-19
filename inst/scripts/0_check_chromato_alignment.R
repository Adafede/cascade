start <- Sys.time()

library(package = baseline, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = docopt, quietly = TRUE)
library(package = future, quietly = TRUE)
library(package = future.apply, quietly = TRUE)
library(package = microshades, quietly = TRUE)
library(package = MSnbase, quietly = TRUE)
library(package = nucleR, quietly = TRUE)
library(package = parallel, quietly = TRUE)
library(package = patchwork, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = progressr, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = xcms, quietly = TRUE)

source(file = "R/baseline_chromatogram.R")
source(file = "R/baseline_chromatograms_progress.R")
source(file = "R/change_intensity_name.R")
source(file = "R/colors.R")
# source(file = "R/dirty_paths_params.R")
source(file = "R/get_params.R")
source(file = "R/improve_signal.R")
source(file = "R/improve_signals_progress.R")
source(file = "R/log_debug.R")
source(file = "R/make_confident.R")
source(file = "R/normalize_chromatograms_list.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/peaks_progress.R")
source(file = "R/plot_histograms.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_plot.R")
source(file = "R/y_as_na.R")

progressr::handlers(global = TRUE)
progressr::handlers("progress")

step <- "processing"
paths <- parse_yaml_paths()
params <- ""
params <- get_params(step = step)

log_debug(
  "This program compares",
  "profiles of different lengths",
  "of in depth annotated extracts"
)
log_debug("Authors: \n", "AR")
log_debug("Contributors: \n", "...")

#' Dirty generic paths and parameters
source(file = "R/dirty_paths_params.R")
#' Specific paths
TOYSET <- "data/source/mzml/10043"
TOYSET <- "data/source/mzml/UHR"
TOYSET <- "data/source/mzml/bitter"
TIME_MIN <- 0.3
TIME_MAX <- 21.3

files <- list.files(
  path = TOYSET,
  # pattern = "10043_enriched_UHR_Pos",
  pattern = "10043_apolar_Pos",
  full.names = TRUE,
  recursive = TRUE
)

files_2 <- list.files(
  path = TOYSET,
  # pattern = "10043_enriched_UHR_Neg",
  pattern = "10043_apolar_Neg",
  full.names = TRUE,
  recursive = TRUE
)

names <- list.files(
  path = TOYSET,
  pattern = "Pos.mzML",
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML",
    replacement = "",
    fixed = TRUE
  )

names_2 <- list.files(
  path = TOYSET,
  pattern = "Neg.mzML",
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML",
    replacement = "",
    fixed = TRUE
  )

# dda_data <- MSnbase::readMSData(
#   files = files,
#   mode = "onDisk",
#   msLevel. = 1
# )

objects <- lapply(files, mzR::openMSfile)
objects_2 <- lapply(files_2, mzR::openMSfile)

chromatograms <- mzR::chromatograms(object = objects[[
  1
]])
chromatograms_2 <- lapply(objects_2, mzR::chromatograms)

chromatograms_all <- chromatograms
chromatograms_all_2 <- purrr::flatten(chromatograms_2)

chromatograms_bpi <- chromatograms_all[c(TRUE, FALSE, FALSE)]

chromatograms_bpi_neg <- chromatograms_all_2[c(TRUE, FALSE, FALSE)]

chromatograms_cad <- chromatograms_all[c(FALSE, FALSE, TRUE)]

chromatograms_pda <- chromatograms_all[c(FALSE, TRUE, FALSE)]

chromatograms_bpi_ready <-
  lapply(chromatograms_bpi, change_intensity_name, name = "BasePeak_0")

chromatograms_bpi_neg_ready <-
  lapply(chromatograms_bpi_neg, change_intensity_name, name = "BasePeak_0")

chromatograms_cad_ready <-
  lapply(chromatograms_cad, change_intensity_name, name = "UV.1_CAD_1_0")

chromatograms_pda_ready <-
  lapply(
    chromatograms_pda,
    change_intensity_name,
    name = "PDA.1_TotalAbsorbance_0"
  )

chromatograms_bpi_improved <-
  improve_signals_progress(chromatograms_bpi_ready)
chromatograms_bpi_neg_improved <-
  improve_signals_progress(chromatograms_bpi_neg_ready)
chromatograms_cad_improved <-
  improve_signals_progress(chromatograms_cad_ready)
chromatograms_pda_improved <-
  improve_signals_progress(chromatograms_pda_ready)

# names(chromatograms_bpi_improved) <- names
# names(chromatograms_bpi_neg_improved) <- names_2
# names(chromatograms_cad_improved) <- names
# names(chromatograms_pda_improved) <- names

bpis_improved <- chromatograms_bpi_improved |>
  normalize_chromatograms_list()

bpis_neg_improved <- chromatograms_bpi_neg_improved |>
  normalize_chromatograms_list()

cads_improved <- chromatograms_cad_improved |>
  normalize_chromatograms_list(shift = CAD_SHIFT)

pdas_improved <- chromatograms_pda_improved |>
  normalize_chromatograms_list(shift = PDA_SHIFT)

chromatograms_bpi_baselined <-
  baseline_chromatograms_progress(chromatograms_bpi_improved)
chromatograms_bpi_neg_baselined <-
  baseline_chromatograms_progress(chromatograms_bpi_neg_improved)
chromatograms_cad_baselined <-
  baseline_chromatograms_progress(chromatograms_cad_improved)
chromatograms_pda_baselined <-
  baseline_chromatograms_progress(chromatograms_pda_improved)

bpis_baselined <- chromatograms_bpi_baselined |>
  normalize_chromatograms_list()
bpis_neg_baselined <- chromatograms_bpi_neg_baselined |>
  normalize_chromatograms_list()
cads_baselined <- chromatograms_cad_baselined |>
  normalize_chromatograms_list()
pdas_baselined <- chromatograms_pda_baselined |>
  normalize_chromatograms_list()

new <- plotly::plot_ly() |>
  plotly::add_lines(
    data = chromatograms_cad_improved[[1]] |>
      normalize_chromatograms_list(shift = CAD_SHIFT, time = TRUE),
    x = ~time,
    y = ~intensity,
    name = "<b> CAD </b>",
    line = list(
      width = 1,
      # dash = "dot",
      color = "e31a1c"
    )
  ) |>
  plotly::add_lines(
    data = chromatograms_pda_improved[[1]] |>
      normalize_chromatograms_list(shift = PDA_SHIFT, time = TRUE),
    x = ~time,
    y = ~intensity,
    name = "<b> PDA </b>",
    line = list(
      width = 1,
      # dash = "dot",
      color = "b2df8a"
    )
  ) |>
  plotly::add_lines(
    data = chromatograms_bpi_improved[[1]] |>
      normalize_chromatograms_list(time = TRUE),
    x = ~time,
    y = ~ -intensity,
    name = "<b> MS Pos </b>",
    line = list(
      width = 1,
      # dash = "dot",
      color = "a6cee3"
    )
  ) |>
  plotly::add_lines(
    data = chromatograms_bpi_neg_improved[[1]] |>
      normalize_chromatograms_list(time = TRUE),
    x = ~time,
    y = ~ -intensity,
    name = "<b> MS Neg </b>",
    line = list(
      width = 1,
      # dash = "dot",
      color = "1f78b4"
    )
  ) |>
  plotly::layout(
    xaxis = list(title = "<b> Time [minutes] </b>"),
    yaxis = list(title = "<b> Normalized Intensity </b>")
  )
new
