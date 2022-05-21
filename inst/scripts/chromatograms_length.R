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

source(file = "R/colors.R")
source(file = "R/get_params.R")
source(file = "R/improve_signal.R")
source(file = "R/improve_signals_progress.R")
source(file = "R/log_debug.R")
source(file = "R/make_confident.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/peaks_progress.R")
source(file = "R/plot_histograms.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_plot.R")
source(file = "R/y_as_na.R")

future::plan(strategy = future::multisession)
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

# dda_data <- MSnbase::readMSData(
#   files = files,
#   mode = "onDisk",
#   msLevel. = 1
# )

objects <- lapply(files, mzR::openMSfile)

chromatograms <- lapply(objects, mzR::chromatograms)

chromatograms_all <- purrr::flatten(chromatograms)

chromatograms_cad <- chromatograms_all[c(FALSE, FALSE, TRUE)]

chromatograms_cad_ready <- lapply(chromatograms_cad, function(x) {
  x |>
    dplyr::select(time,
      intensity = UV.1_CAD_1_0
    )
})

chromatograms_cad_improved <-
  improve_signals_progress(chromatograms_cad_ready)

names(chromatograms_cad_improved) <- names

cads_improved <-
  dplyr::bind_rows(chromatograms_cad_improved, .id = "id") |>
  dplyr::mutate(time = time + CAD_SHIFT) |>
  dplyr::group_by(id) |>
  dplyr::mutate(time_2 = max(time)) |>
  dplyr::mutate(time = time / time_2)

cad_plot <- plotly::plot_ly(
  data = cads_improved,
  x = ~time,
  y = ~intensity,
  color = ~id,
  colors = "Spectral",
  type = "scatter",
  mode = "lines",
  line = list(width = 0.5),
  legendgroup = ~id
) |>
  plotly::layout(
    annotations = list(
      x = 0.95,
      y = 0.95,
      xref = "paper",
      yref = "paper",
      text = "CAD",
      showarrow = FALSE
    )
  )

chromatograms_cad_baselined <- chromatograms_cad_improved

for (i in seq_along(seq_len(length(chromatograms_cad_improved)))) {
  intensity <- chromatograms_cad_improved[[i]]$intensity

  intensity[is.na(intensity)] <- 0

  intensity_baseline <- baseline(
    spectra = t(intensity),
    method = "peakDetection"
  )

  intensity_new <- t(intensity_baseline@corrected) |>
    data.table::data.table()

  chromatograms_cad_baselined[[i]]$intensity <- intensity_new$V1
}

cads_baselined <-
  dplyr::bind_rows(chromatograms_cad_baselined, .id = "id") |>
  dplyr::mutate(intensity = intensity / max(intensity)) |>
  dplyr::group_by(id) |>
  dplyr::mutate(time_2 = max(time)) |>
  dplyr::mutate(time = time / time_2)

new_new <- plotly::plot_ly(
  data = cads_baselined,
  # |> dplyr::filter(grepl(pattern = "M", x = id)),
  x = ~time,
  y = ~intensity,
  color = ~id,
  colors = "Spectral",
  type = "scatter",
  mode = "lines",
  line = list(width = 1),
  legendgroup = ~id
) |>
  plotly::layout(
    annotations = list(
      x = 0.95,
      y = 0.95,
      xref = "paper",
      yref = "paper",
      text = "CAD",
      showarrow = FALSE
    )
  )

uhr <- plotly::plot_ly() |>
  plotly::add_lines(
    data = cads_baselined |> dplyr::filter(grepl(pattern = "UHR", x = id)),
    x = ~time,
    y = ~intensity,
    name = "very long"
  )
semi <- plotly::plot_ly() |>
  plotly::add_lines(
    data = cads_baselined |> dplyr::filter(grepl(pattern = "21", x = id)),
    x = ~time,
    y = ~intensity,
    name = "long"
  )
short <- plotly::plot_ly() |>
  plotly::add_lines(
    data = cads_baselined |> dplyr::filter(grepl(pattern = "09_P", x = id)),
    x = ~time,
    y = ~intensity,
    name = "short"
  )

new_new_new <- plotly::subplot(short,
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
    yaxis2 = list(range = c(0, 0.2)),
    yaxis3 = list(range = c(0, 0.2))
  )
new_new_new

#' WIP

end <- Sys.time()

log_debug("Script finished in", format(end - start))
