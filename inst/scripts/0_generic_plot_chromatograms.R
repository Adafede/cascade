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
source(file = "R/get_params.R")
source(file = "R/improve_signal.R")
source(file = "R/improve_signals_progress.R")
source(file = "R/log_debug.R")
source(file = "R/make_confident.R")
source(file = "R/normalize_chromatograms_list.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
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

#' Paths
TOYSET <- "~/data/qcmix/"

#' Generic parameters
WORKERS <- params$workers

#' Parameters for LC alignment
TIME_MIN <- params$chromato$time$min
TIME_MAX <- params$chromato$time$max
CAD_SHIFT <- params$chromato$shift$cad
PDA_SHIFT <- params$chromato$shift$pda
ESTIMATED_SOLUBLITIY_LIMIT <- params$misc$solubility$limit

#' Parameters for signal improvement
FOURRIER_COMPONENTS <- params$signal$fourrier$components
FREQUENCY <- params$signal$frequency
RESAMPLE <- params$signal$resample

#' Parameters adapted from Excel sheet from paper shortDOI: 10/ghmvhz
sigma <- params$signal$sigma
k2 <- sigma / params$signal$k2 # 30
k4 <- sigma / params$signal$k4 # 200
smoothing_width <- params$signal$smoothing
baseline_adjust <- params$signal$baseline

#' Parameters related to MS/CAD
INTENSITY <- params$chromato$intensity$ms1$min
PEAK_SIMILARITY <- params$chromato$peak$similarity$filter
PEAK_SIMILARITY_PREFILTER <-
  params$chromato$peak$similarity$prefilter
RT_TOL <- params$chromato$peak$tolerance$rt
PPM <- params$chromato$peak$tolerance$ppm

#' Parameters for annotation
CONFIDENCE_SCORE_MIN <- params$annotation$confidence$min

files <- list.files(
  path = TOYSET,
  pattern = "Pos.mzML.gz",
  full.names = TRUE,
  recursive = TRUE
)

files_2 <- list.files(
  path = TOYSET,
  pattern = "Neg.mzML.gz",
  full.names = TRUE,
  recursive = TRUE
)

names <- list.files(
  path = TOYSET,
  pattern = "Pos.mzML.gz",
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML.gz",
    replacement = "",
    fixed = TRUE
  )

names_2 <- list.files(
  path = TOYSET,
  pattern = "Neg.mzML.gz",
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML.gz",
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

chromatograms <- lapply(objects, mzR::chromatograms)
chromatograms_2 <- lapply(objects_2, mzR::chromatograms)

chromatograms_all <- purrr::flatten(chromatograms)
chromatograms_all_2 <- purrr::flatten(chromatograms_2)


chromatograms_bpi <- chromatograms_all[c(TRUE, FALSE, FALSE)]

chromatograms_bpi_neg <- chromatograms_all_2[c(TRUE, FALSE, FALSE)]

chromatograms_cad <- chromatograms_all[c(FALSE, FALSE, TRUE)]

chromatograms_pda <- chromatograms_all[c(FALSE, TRUE, FALSE)]

chromatograms_bpi_ready <-
  lapply(chromatograms_bpi, change_intensity_name, "BasePeak_0")

chromatograms_bpi_neg_ready <-
  lapply(chromatograms_bpi_neg, change_intensity_name, "BasePeak_0")

chromatograms_cad_ready <-
  lapply(chromatograms_cad, change_intensity_name, "UV.1_CAD_1_0")

chromatograms_pda_ready <-
  lapply(
    chromatograms_pda,
    change_intensity_name,
    "PDA.1_TotalAbsorbance_0"
  )

chromatograms_bpi_improved <-
  improve_signals_progress(chromatograms_bpi_ready)
chromatograms_bpi_neg_improved <-
  improve_signals_progress(chromatograms_bpi_neg_ready)
chromatograms_cad_improved <-
  improve_signals_progress(chromatograms_cad_ready)
chromatograms_pda_improved <-
  improve_signals_progress(chromatograms_pda_ready)

names(chromatograms_bpi_improved) <- names
names(chromatograms_bpi_neg_improved) <- names_2
names(chromatograms_cad_improved) <- names
names(chromatograms_pda_improved) <- names

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

uhr <- plotly::plot_ly() |>
  plotly::add_lines(
    data = cads_baselined |>
      dplyr::filter(grepl(pattern = "UHR", x = id)),
    x = ~time,
    y = ~intensity,
    name = "very long",
    line = list(
      width = 1,
      color = "ff7f00"
    )
  )
semi <- plotly::plot_ly() |>
  plotly::add_lines(
    data = cads_baselined |>
      dplyr::filter(grepl(pattern = "21", x = id)),
    x = ~time,
    y = ~intensity,
    name = "long",
    line = list(
      width = 1,
      color = "cab2d6"
    )
  )
short <- plotly::plot_ly() |>
  plotly::add_lines(
    data = cads_baselined |>
      dplyr::filter(grepl(pattern = "09_P", x = id)),
    x = ~time,
    y = ~intensity,
    name = "short",
    line = list(
      width = 1,
      color = "b15928"
    )
  )

new_new_new <- plotly::subplot(
  short,
  semi,
  uhr,
  shareX = TRUE,
  shareY = TRUE,
  titleY = FALSE,
  titleX = FALSE,
  nrows = 3
) |>
  plotly::layout(
    yaxis = list(range = c(0, 0.2)),
    yaxis2 = list(
      title = "Intensity",
      range = c(0, 0.2)
    ),
    yaxis3 = list(range = c(0, 0.2)),
    xaxis = list(title = "Time")
  )
new_new_new

new <- plotly::plot_ly() |>
  plotly::add_lines(
    data = cads_baselined |>
      dplyr::filter(grepl(pattern = "UHR", x = id)),
    x = ~time,
    y = ~intensity,
    name = "CAD",
    line = list(
      width = 1,
      # dash = "dot",
      color = "e31a1c"
    )
  ) |>
  plotly::add_lines(
    data = pdas_baselined |>
      dplyr::filter(grepl(pattern = "UHR", x = id)),
    x = ~time,
    y = ~intensity,
    name = "PDA",
    line = list(
      width = 1,
      # dash = "dot",
      color = "b2df8a"
    )
  ) |>
  plotly::add_lines(
    data = bpis_baselined |>
      dplyr::filter(grepl(pattern = "UHR", x = id)),
    x = ~time,
    y = ~ -intensity,
    name = "MS Pos",
    line = list(
      width = 1,
      # dash = "dot",
      color = "a6cee3"
    )
  ) |>
  plotly::add_lines(
    data = bpis_neg_baselined |>
      dplyr::filter(grepl(pattern = "UHR", x = id)),
    x = ~time,
    y = ~ -intensity,
    name = "MS Neg",
    line = list(
      width = 1,
      # dash = "dot",
      color = "1f78b4"
    )
  ) |>
  plotly::layout(
    xaxis = list(title = "Time"),
    yaxis = list(title = "Intensity")
  )
new

new_zoom <- new |>
  plotly::layout(xaxis = list(range = c(0.18, 0.33)),
                 yaxis = list(range = c(-0.3, 0.3)))
new_zoom

#' WIP

end <- Sys.time()

log_debug("Script finished in", format(end - start))
