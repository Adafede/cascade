start <- Sys.time()

library(package = baseline, quietly = TRUE)
library(package = data.table, quietly = TRUE)
library(package = docopt, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = future, quietly = TRUE)
library(package = future.apply, quietly = TRUE)
library(package = MSnbase, quietly = TRUE)
library(package = nucleR, quietly = TRUE)
library(package = parallel, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = pracma, quietly = TRUE)
library(package = progressr, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = xcms, quietly = TRUE)

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

progressr::handlers(global = TRUE)
progressr::handlers("progress")

step <- "processing"
paths <- parse_yaml_paths()
params <- ""
params <- get_params(step = step)

log_debug(
  "This program performs",
  "Quantitative and Qualitative Contextualization",
  "of in depth annotated extracts"
)
log_debug("Authors: \n", "AR")
log_debug("Contributors: \n", "...")

#' Dirty generic paths and parameters
source(file = "R/dirty_paths_params.R")
#' Specific paths

log_debug(x = "listing files")
files <- list.files(
  path = TOYSET,
  pattern = paste0(params$filename$mzml, ".mzML"),
  full.names = TRUE,
  recursive = TRUE
)

# files <- files[grepl(pattern = "M_17|M_28|M_36|M_40|M_47|M_57|M_67", x = files)]

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

# names <- names[grepl(pattern = "M_17|M_28|M_36|M_40|M_47|M_57|M_67", x = names)]

log_debug(x = "loading feature table")
# feature_table <- read_features(id = GNPS_JOB)
feature_table <- readr::read_delim(file = FEATURES)

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

log_debug(x = "preparing feature list ...")
df_features <- feature_table |>
  prepare_features()

if (params$signal$detector$bpi == TRUE) {
  chromatograms_list_bpi <- preprocess_chromatograms(
    detector = "bpi",
    list = chromatograms_all[c(TRUE, FALSE, FALSE)],
    signal_name = "BasePeak_0",
    shift = 0
  )
  peaks_prelist_bpi <- preprocess_peaks(
    detector = "bpi",
    list = chromatograms_list_bpi$chromatograms_improved,
    df_long = chromatograms_list_bpi$chromatograms_improved_long
  )
  detector <- "bpi"
  peaks_list_bpi <- process_peaks(detector = "bpi")
}
if (params$signal$detector$cad == TRUE) {
  chromatograms_list_cad <- preprocess_chromatograms()
  peaks_prelist_cad <- preprocess_peaks()
  detector <- "cad"
  peaks_list_cad <- process_peaks()
}
if (params$signal$detector$pda == TRUE) {
  #' TODO check if change for PDA doing baselining on normal and not improved
  chromatograms_list_pda <-
    preprocess_chromatograms(
      detector = "pda",
      list = chromatograms_all[c(FALSE, TRUE, FALSE)],
      signal_name = "PDA.1_TotalAbsorbance_0",
      shift = PDA_SHIFT
    )
  peaks_prelist_pda <- preprocess_peaks(
    detector = "pda",
    list = chromatograms_list_pda$chromatograms_improved,
    df_long = chromatograms_list_pda$chromatograms_improved_long
  )
  detector <- "pda"
  peaks_list_pda <- process_peaks(detector = "pda")
}

#' TODO ADD PEAK PICKING COMPARISON
detector <- "cad"

peaks_original <-
  preprocess_peaks(
    list = chromatograms_list_cad$chromatograms_original,
    df_long = chromatograms_list_cad$chromatograms_original_long
  )
peaks_improved <-
  preprocess_peaks(
    list = chromatograms_list_cad$chromatograms_improved,
    df_long = chromatograms_list_cad$chromatograms_improved_long
  )

#' WIP
suite_1_1 <- chromatograms_list_cad$chromatograms_original_long |>
  dplyr::bind_rows() |>
  # dplyr::filter(grepl(pattern = "V_03_2_01", x = id)) |>
  dplyr::filter(row_number() %% 10 == 1)

suite_1_2 <- peaks_original$list_df_features_with_peaks_long |>
  dplyr::bind_rows() # |>
# dplyr::filter(grepl(pattern = "V_03_2_01", x = id))

suite_2_1 <- chromatograms_list_cad$chromatograms_improved_long |>
  dplyr::bind_rows() # |>
# dplyr::filter(grepl(pattern = "V_03_2_01", x = id))

suite_2_2 <- peaks_improved$list_df_features_with_peaks_long |>
  dplyr::bind_rows() # |>
# dplyr::filter(grepl(pattern = "V_03_2_01", x = id))

suite_3_1 <- chromatograms_list_cad$chromatograms_baselined_long |>
  dplyr::bind_rows() # |>
# dplyr::filter(grepl(pattern = "V_03_2_01", x = id))

suite_3_2 <- peaks_prelist_cad$list_df_features_with_peaks_long |>
  dplyr::bind_rows() # |>
# dplyr::filter(grepl(pattern = "V_03_2_01", x = id))


f_1 <- approxfun(
  x = chromatograms_list_cad$chromatograms_original_long |>
    dplyr::bind_rows() |>
    # dplyr::filter(grepl(pattern = "V_03_2_01", x = id)) |>
    dplyr::pull(time),
  y = chromatograms_list_cad$chromatograms_original_long |>
    dplyr::bind_rows() |>
    # dplyr::filter(grepl(pattern = "V_03_2_01", x = id)) |>
    dplyr::pull(intensity)
)

detection_before <- plotly::plot_ly(suite_1_1) |>
  plotly::add_trace(
    data = suite_1_1,
    x =  ~time,
    y =  ~intensity,
    type = "scatter",
    mode = "line",
    name = "signal",
    line = list(color = "1f78b4")
  ) |>
  # plotly::add_trace(
  #   data = suite_1_2,
  #   x = ~rt_apex,
  #   y = ~peak_max,
  #   yaxis = "y2",
  #   type = "scatter",
  #   marker = list(
  #     color = "ff7f00",
  #     symbol = "star"
  #   ),
  #   name = "detected maximum",
  #   line = list(color = "1f78b4", width = 0)
  # ) |>
  # plotly::add_trace(
  #   data = suite_1_2,
  #   x = ~rt_min,
  #   y = ~ f_1(rt_min),
  #   yaxis = "y2",
  #   type = "scatter",
  #   marker = list(
  #     color = "33a02c",
  #     symbol = "triangle-right"
  #   ),
  #   name = "detected minimum (start)",
  #   line = list(color = "1f78b4", width = 0)
  # ) |>
  # plotly::add_trace(
  #   data = suite_1_2,
  #   x = ~rt_max,
  #   y = ~ f_1(rt_max),
  #   yaxis = "y2",
  #   type = "scatter",
  #   marker = list(
  #     color = "33a02c",
  #     symbol = "triangle-left"
  #   ),
  #   name = "detected minimum (end)",
  #   line = list(color = "1f78b4", width = 0)
  # ) |>
  plotly::layout(
    yaxis = list(
      # range = c(0, 1),
      title = "Normalized Intensity",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    ),
    yaxis2 = list(
      # range = c(0, 1),
      title = "",
      showticklabels = FALSE,
      overlaying = "y",
      side = "right"
    ),
    xaxis = list(
      range = c(0.5, 127),
      title = "Time",
      showticklabels = FALSE,
      overlaying = "y",
      side = "right"
    ),
    showlegend = FALSE
  )

detection_before

f_2 <- approxfun(
  x = chromatograms_list_cad$chromatograms_improved_long |>
    dplyr::bind_rows() |>
    # dplyr::filter(grepl(pattern = "V_03_2_01", x = id)) |>
    dplyr::pull(time),
  y = chromatograms_list_cad$chromatograms_improved_long |>
    dplyr::bind_rows() |>
    # dplyr::filter(grepl(pattern = "V_03_2_01", x = id)) |>
    dplyr::pull(intensity)
)

detection_after <- plotly::plot_ly(suite_2_1) |>
  plotly::add_trace(
    suite_2_1,
    x =  ~time,
    y =  ~intensity,
    type = "scatter",
    mode = "line",
    name = "signal",
    line = list(color = "1f78b4")
  ) |>
  plotly::add_trace(
    data = suite_2_2,
    x = ~rt_apex,
    y = ~peak_max,
    yaxis = "y2",
    type = "scatter",
    marker = list(
      color = "ff7f00",
      symbol = "star"
    ),
    name = "detected maximum",
    line = list(color = "1f78b4", width = 0)
  ) |>
  plotly::add_trace(
    data = suite_2_2,
    x = ~rt_min,
    y = ~ f_2(rt_min),
    yaxis = "y2",
    type = "scatter",
    marker = list(
      color = "33a02c",
      symbol = "triangle-right"
    ),
    name = "detected minimum (start)",
    line = list(color = "1f78b4", width = 0)
  ) |>
  plotly::add_trace(
    data = suite_2_2,
    x = ~rt_max,
    y = ~ f_2(rt_max),
    yaxis = "y2",
    type = "scatter",
    marker = list(
      color = "33a02c",
      symbol = "triangle-left"
    ),
    name = "detected minimum (end)",
    line = list(color = "1f78b4", width = 0)
  ) |>
  plotly::layout(
    yaxis = list(
      range = c(0, 1),
      title = "Normalized Intensity",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    ),
    yaxis2 = list(
      range = c(0, 1),
      title = "",
      showticklabels = FALSE,
      overlaying = "y",
      side = "right"
    ),
    xaxis = list(
      range = c(0.5, 127),
      title = "Time",
      showticklabels = FALSE,
      overlaying = "y",
      side = "right"
    )
    # showlegend = FALSE
  )

detection_after

# f_3 <- approxfun(
#   x = chromatograms_list_cad$chromatograms_baselined_long |>
#     dplyr::bind_rows() |>
#     # dplyr::filter(grepl(pattern = "V_03_2_01", x = id)) |>
#     dplyr::pull(time),
#   y = chromatograms_list_cad$chromatograms_baselined_long |>
#     dplyr::bind_rows() |>
#     # dplyr::filter(grepl(pattern = "V_03_2_01", x = id)) |>
#     dplyr::pull(intensity)
# )
#
# plotly::plot_ly(
#   suite_3_1,
#   x =  ~ time,
#   y =  ~ intensity,
#   type = "scatter",
#   mode = "line"
# ) |>
#   plotly::add_trace(
#     data = suite_3_2,
#     x = ~ rt_apex,
#     y = ~ peak_max,
#     yaxis = "y2",
#     type = "scatter",
#     marker = list(color = "ff7f00")
#   ) |>
#   plotly::add_trace(
#     data = suite_3_2,
#     x = ~ rt_min,
#     y = ~ f_3(rt_min),
#     yaxis = "y2",
#     type = "scatter",
#     marker = list(color = "green")
#   ) |>
#   plotly::add_trace(
#     data = suite_3_2,
#     x = ~ rt_max,
#     y = ~ f_3(rt_max),
#     yaxis = "y2",
#     type = "scatter",
#     marker = list(color = "green")
#   ) |>
#   plotly::layout(
#     yaxis = list(range = c(0, 1)),
#     yaxis2 = list(
#       range = c(0, 1),
#       overlaying = "y",
#       side = "right"
#     )
#   )

# Sa_2 <- chromatograms_list_cad$chromatograms_baselined_long |>
#   dplyr::bind_rows() |>
#   tibble::column_to_rownames(var = "time") |>
#   dplyr::select("666" = intensity) |>
#   as.matrix()
#
# new.ts_2 <- rownames(Sa_2) |> as.numeric()
# new.lambdas_2 <- colnames(Sa_2) |> as.numeric()
#
# dat.pr <- list(Sa_2)
# names(dat.pr) <- "666"
#
# pks_gau <-
#   get_peaks(
#     chrom_list = dat.pr,
#     lambdas = new.lambdas_2,
#     sd.max = 25,
#     max.iter = 10000,
#     fit = "gauss"
#   )
# plot(pks_gau, index = 1, lambda = '666')

# plotly::save_image(
#   p = detection_before,
#   file = "data/chromatograms/detection_before.pdf",
#   height = 450,
#   width = 800,
#   scale = 3
# )
# plotly::save_image(
#   p = detection_after,
#   file = "data/chromatograms/detection_after.pdf",
#   height = 450,
#   width = 800,
#   scale = 3
# )

end <- Sys.time()

log_debug("Script finished in", format(end - start))
