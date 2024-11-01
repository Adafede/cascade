start <- Sys.time()

pkgload::load_all()

message(
  "This program performs",
  "Quantitative and Qualitative Contextualization",
  "of in depth annotated extracts"
)
message("Authors: \n", "AR")
message("Contributors: \n", "...")

#' Specific paths
AREA_MIN <- 0.005
CAD_SHIFT <- 0.05
INTENSITY_MS_MIN <- 10000
TIME_MIN <- 0.7
TIME_MAX <- 35.2
FILE_POSITIVE <- "data/source/mzml/210619_AR_06_V_03_2_01.mzML"
FEATURES <- "~/Documents/papers/sapid/sapere_tmp/extract_mzmine/extract.csv"

name <- FILE_POSITIVE |>
  gsub(pattern = ".*/", replacement = "") |>
  gsub(pattern = "[0-9]{8}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML",
    replacement = "",
    fixed = TRUE
  )

message("loading feature table")
feature_table <- readr::read_delim(file = FEATURES)

message("opening raw files objects and extracting chromatograms")
chromatograms_all <- load_chromatograms(show_example = TRUE)

message("preparing feature list ...")
df_features <- feature_table |>
  prepare_features(min_intensity = INTENSITY_MS_MIN, name = name)

chromatograms_list_cad <- preprocess_chromatograms(name = name)
peaks_prelist_cad <- preprocess_peaks(area_min = AREA_MIN, name = name)

detector <- "cad"

peaks_original <-
  preprocess_peaks(
    list = chromatograms_list_cad$chromatograms_original,
    df_long = chromatograms_list_cad$chromatograms_original_long |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    area_min = AREA_MIN,
    name = name
  )
peaks_improved <-
  preprocess_peaks(
    list = chromatograms_list_cad$chromatograms_improved,
    df_long = chromatograms_list_cad$chromatograms_improved_long |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    area_min = AREA_MIN,
    name = name
  )
peaks_baselined <-
  preprocess_peaks(
    list = chromatograms_list_cad$chromatograms_baselined,
    df_long = chromatograms_list_cad$chromatograms_baselined_long |>
      tidytable::mutate(intensity = intensity / max(intensity)),
    area_min = AREA_MIN, name = name
  )

suite_1_1 <- chromatograms_list_cad$chromatograms_original_long |>
  tidytable::bind_rows() |>
  tidytable::mutate(intensity = intensity / max(intensity)) |>
  tidytable::filter(row_number() %% 10 == 1)

suite_1_2 <- peaks_original$list_df_features_with_peaks_long |>
  tidytable::bind_rows()

suite_2_1 <- chromatograms_list_cad$chromatograms_improved_long |>
  tidytable::bind_rows() |>
  tidytable::mutate(intensity = intensity / max(intensity))

suite_2_2 <- peaks_improved$list_df_features_with_peaks_long |>
  tidytable::bind_rows() |>
  tidytable::mutate(
    intensity = intensity_max / max(intensity_max),
    peak_max = peak_max / max(peak_max)
  )

suite_3_1 <- chromatograms_list_cad$chromatograms_baselined_long |>
  tidytable::bind_rows()

suite_3_2 <- peaks_prelist_cad$list_df_features_with_peaks_long |>
  tidytable::bind_rows() |>
  tidytable::mutate(
    intensity = intensity_max / max(intensity_max),
    peak_max = peak_max / max(peak_max)
  )


f_1 <- approxfun(
  x = chromatograms_list_cad$chromatograms_original_long |>
    tidytable::bind_rows() |>
    tidytable::mutate(intensity = intensity / max(intensity)) |>
    tidytable::pull(time),
  y = chromatograms_list_cad$chromatograms_original_long |>
    tidytable::bind_rows() |>
    tidytable::mutate(intensity = intensity / max(intensity)) |>
    tidytable::pull(intensity)
)

detection_before <- plotly::plot_ly(suite_1_1) |>
  plotly::add_trace(
    data = suite_1_1,
    x =  ~time,
    y =  ~intensity,
    type = "scatter",
    mode = "line",
    name = "signal",
    line = list(color = "1f78b4", width = 1)
  ) |>
  plotly::add_trace(
    data = suite_1_2,
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
    data = suite_1_2,
    x = ~rt_min,
    y = ~ f_1(rt_min),
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
    data = suite_1_2,
    x = ~rt_max,
    y = ~ f_1(rt_max),
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
      title = "Normalized Intensity",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    ),
    yaxis2 = list(
      title = "",
      showticklabels = FALSE,
      overlaying = "y",
      side = "right"
    ),
    xaxis = list(
      title = "Time [min]",
      showticklabels = FALSE,
      showgrid = FALSE,
      overlaying = "y",
      side = "right"
    ),
    showlegend = FALSE
  )

detection_before

f_2 <- approxfun(
  x = chromatograms_list_cad$chromatograms_baselined_long |>
    tidytable::bind_rows() |>
    tidytable::mutate(intensity = intensity / max(intensity)) |>
    tidytable::pull(time),
  y = chromatograms_list_cad$chromatograms_baselined_long |>
    tidytable::bind_rows() |>
    tidytable::mutate(intensity = intensity / max(intensity)) |>
    tidytable::pull(intensity)
)

detection_after <- suite_3_1 |>
  plotly::plot_ly() |>
  plotly::add_trace(
    suite_3_1,
    x =  ~time,
    y =  ~intensity,
    type = "scatter",
    mode = "line",
    name = "signal",
    line = list(color = "1f78b4", width = 1)
  ) |>
  plotly::add_trace(
    data = suite_3_2,
    x = ~rt_apex,
    y = ~peak_max,
    yaxis = "y2",
    type = "scatter",
    marker = list(color = "ff7f00", symbol = "star"),
    name = "detected maximum",
    line = list(color = "1f78b4", width = 0)
  ) |>
  plotly::add_trace(
    data = suite_3_2,
    x = ~rt_min,
    y = ~ f_2(rt_min),
    yaxis = "y2",
    type = "scatter",
    marker = list(color = "33a02c", symbol = "triangle-right"),
    name = "detected minimum (start)",
    line = list(color = "1f78b4", width = 0)
  ) |>
  plotly::add_trace(
    data = suite_3_2,
    x = ~rt_max,
    y = ~ f_2(rt_max),
    yaxis = "y2",
    type = "scatter",
    marker = list(color = "33a02c", symbol = "triangle-left"),
    name = "detected minimum (end)",
    line = list(color = "1f78b4", width = 0)
  ) |>
  plotly::layout(
    yaxis = list(
      title = "Normalized Intensity",
      zeroline = TRUE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    ),
    yaxis2 = list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE,
      overlaying = "y",
      side = "right"
    ),
    xaxis = list(
      title = "Time [min]",
      showticklabels = FALSE,
      showgrid = FALSE,
      overlaying = "y",
      side = "right"
    ),
    showlegend = FALSE
  )
detection_after

end <- Sys.time()

message("Script finished in ", format(end - start))
