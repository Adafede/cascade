start <- Sys.time()

library(package = baseline, quietly = TRUE)
library(package = data.table, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = docopt, quietly = TRUE)
library(package = microshades, quietly = TRUE)
library(package = MSnbase, quietly = TRUE)
library(package = nucleR, quietly = TRUE)
library(package = patchwork, quietly = TRUE)
library(package = parallel, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = pracma, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = xcms, quietly = TRUE)

source(file = "R/colors.R")
source(file = "R/compare_peaks.R")
source(file = "R/extract_peak.R")
source(file = "R/get_gnps.R")
source(file = "R/get_params.R")
source(file = "R/improve_signal.R")
source(file = "R/log_debug.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/plot_histograms.R")
source(file = "R/prepare_plot.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_hierarchy_preparation.R")
source(file = "R/y_as_na.R")

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

#' Paths
ANNOTATIONS <-
  "~/git/tima-r/inst/extdata/processed/211230_110947/20211227_10043.tsv.gz"
GNPS_JOB <- "97d7c50031a84b9ba2747e883a5907cd"
TOYSET <- "~/data/20210701_10043/fractions"

#' Generic parameters
WORKERS <- 10

#' Parameters for LC alignment
TIME_MIN <- 0 ## minutes
TIME_MAX <- 32 ## minutes
CAD_SHIFT <- 0.055 ## minutes
PDA_SHIFT <- 0.090 ## minutes
ESTIMATED_SOLUBLITIY_LIMIT <- 58

#' Parameters for signal improvement
FOURRIER_COMPONENTS <- 0.05 ## of total
FREQUENCY <- 20 ## Hz
RESAMPLE <- 1 / 10 ## points

#' Parameters adapted from Excel sheet from paper shortDOI: 10/ghmvhz
sigma <- 0.05
k2 <- sigma / 250 # 30
k4 <- sigma / 1250000 # 200
smoothing_width <- 5
baseline_adjust <- 0

#' Parameters related to MS/CAD
INTENSITY <- 1E5
PEAK_SIMILARITY <- 0.9
PEAK_SIMILARITY_PREFILTER <- 0.6
PPM <- 10

#' Parameters for annotation
CONFIDENCE_SCORE_MIN <- 0.5

future::plan(strategy = future::multiprocess(workers = WORKERS))

files <- list.files(
  path = TOYSET,
  pattern = ".mzML.gz",
  full.names = TRUE,
  recursive = TRUE
)

names <- list.files(
  path = TOYSET,
  pattern = ".mzML.gz",
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML.gz",
    replacement = "",
    fixed = TRUE
  )

annotations <-
  readr::read_delim(file = ANNOTATIONS)

feature_table <-
  read_features(id = GNPS_JOB)

dda_data <- MSnbase::readMSData(
  files = files,
  mode = "onDisk",
  msLevel. = 1
)

objects <- list()

for (i in seq_along(files)) {
  objects[[i]] <-
    mzR::openMSfile(files[i])
}

chromatograms <- list()

for (i in seq_along(objects)) {
  chromatograms[[i]] <-
    mzR::chromatograms(objects[[i]])
}

chromatograms_all <- purrr::flatten(chromatograms)

# chromatograms_bpi <- chromatograms_all[c(TRUE, FALSE, FALSE)]

# chromatograms_pda <- chromatograms_all[c(FALSE, TRUE, FALSE)]

chromatograms_cad <- chromatograms_all[c(FALSE, FALSE, TRUE)]

# chromatograms_bpi_improved <- list()

# for (i in seq_along(1:length(chromatograms_bpi))) {
#   chromatograms_bpi_improved[[i]] <-
#     improve_signal(df = chromatograms_bpi[[i]] |>
#       dplyr::select(time, intensity = BasePeak_0))
# }

# names(chromatograms_bpi_improved) <- names

# bpis_improved <-
#   dplyr::bind_rows(chromatograms_bpi_improved, .id = "id") |>
#   dplyr::mutate(time = time)

# chromatograms_pda_improved <- list()

# for (i in seq_along(1:length(chromatograms_pda))) {
#   chromatograms_pda_improved[[i]] <-
#     improve_signal(df = chromatograms_pda[[i]] |>
#       dplyr::select(time, intensity = PDA.1_TotalAbsorbance_0))
# }

# names(chromatograms_pda_improved) <- names

# pdas_improved <-
#   dplyr::bind_rows(chromatograms_pda_improved, .id = "id") |>
#   dplyr::mutate(time = time + PDA_SHIFT)

chromatograms_cad_improved <- list()

for (i in seq_along(1:length(chromatograms_cad))) {
  chromatograms_cad_improved[[i]] <-
    improve_signal(df = chromatograms_cad[[i]] |>
      dplyr::select(time, intensity = UV.1_CAD_1_0))
}

names(chromatograms_cad_improved) <- names

cads_improved <-
  dplyr::bind_rows(chromatograms_cad_improved, .id = "id") |>
  dplyr::mutate(time = time + CAD_SHIFT)

# bpi_plot <- plotly::plot_ly(
#   data = bpis_improved,
#   # |> dplyr::filter(grepl(pattern = "M", x = id)),
#   x = ~time,
#   y = ~intensity,
#   color = ~id,
#   colors = "Spectral",
#   type = "scatter",
#   mode = "lines",
#   line = list(width = 0.5),
#   legendgroup = ~id
# ) |>
#   plotly::layout(
#     annotations = list(
#       x = 0.95,
#       y = 0.95,
#       xref = "paper",
#       yref = "paper",
#       text = "MS (POS)",
#       showarrow = FALSE
#     )
#   )

# pda_plot <- plotly::plot_ly(
#   data = pdas_improved,
#   # |> dplyr::filter(grepl(pattern = "M", x = id)),
#   x = ~time,
#   y = ~intensity,
#   color = ~id,
#   colors = "Spectral",
#   type = "scatter",
#   mode = "lines",
#   line = list(width = 0.5),
#   legendgroup = ~id
# ) |>
#   plotly::layout(
#     annotations = list(
#       x = 0.95,
#       y = 0.95,
#       xref = "paper",
#       yref = "paper",
#       text = "UV (200-500nm)",
#       showarrow = FALSE
#     )
#   )

cad_plot <- plotly::plot_ly(
  data = cads_improved,
  # |> dplyr::filter(grepl(pattern = "M", x = id)),
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

# comparison <-
#   plotly::subplot(bpi_plot,
#     pda_plot,
#     cad_plot,
#     nrows = 3,
#     shareX = TRUE
#   )
#
# comparison

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
  dplyr::mutate(intensity = intensity / max(intensity))

new_new <- plotly::plot_ly(
  data = cads_baselined,
  # |> dplyr::filter(grepl(pattern = "M", x = id)),
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

# comparison_2 <-
#   plotly::subplot(bpi_plot,
#     pda_plot,
#     new_new,
#     nrows = 3,
#     shareX = TRUE
#   )
#
# comparison_2

peaks_cad <- list()

for (i in seq_along(seq_len(length(chromatograms_cad_baselined)))) {
  # plot(cads_baselined$intensity, type = "l", col = "navy")
  # grid()
  x <-
    pracma::findpeaks(
      chromatograms_cad_baselined[[i]]$intensity,
      npeaks = 2000,
      threshold = 0.01,
      sortstr = TRUE
    )
  # points(x[, 2], x[, 1], pch = 20, col = "maroon") ## End(Not run)

  peaks <- data.frame(x) |>
    dplyr::mutate(
      peak_id = dplyr::row_number(),
      peak_max = X1,
      rt_apex = chromatograms_cad_baselined[[i]]$time[X2],
      rt_min = chromatograms_cad_baselined[[i]]$time[X3],
      rt_max = chromatograms_cad_baselined[[i]]$time[X4]
    ) |>
    data.table::data.table()

  peaks_cad[[i]] <- peaks
}

names(peaks_cad) <- names

peaks_all <- dplyr::bind_rows(peaks_cad, .id = "id")

cads_baselined <- cads_baselined |>
  dplyr::mutate(rt_1 = time, rt_2 = time) |>
  data.table::data.table()

cat("setting joining keys \n")
data.table::setkey(peaks_all, rt_min, rt_max)
data.table::setkey(cads_baselined, rt_1, rt_2)

cat("joining within given rt tolerance \n")
df <- data.table::foverlaps(peaks_all, cads_baselined) |>
  dplyr::filter(id == i.id) |>
  dplyr::group_by(peak_id, id) |>
  dplyr::mutate(integral = sum(intensity)) |>
  dplyr::ungroup() |>
  dplyr::distinct(peak_id, id, peak_max, rt_apex, rt_min, rt_max, integral) |>
  dplyr::group_by(id) |>
  dplyr::filter(integral >= 0.01 * integral / sum(integral)) |>
  data.table::data.table()

ms1_best_candidate <- annotations |>
  dplyr::mutate(
    smiles_2D = ifelse(
      test = score_final >= CONFIDENCE_SCORE_MIN,
      yes = smiles_2D,
      no = NA
    )
  ) |>
  dplyr::mutate_all(list(~ gsub(
    pattern = "\\|.*",
    replacement = "",
    x = .x
  ))) |>
  splitstackshape::cSplit("best_candidate", sep = "§") |>
  dplyr::distinct(
    feature_id,
    mz,
    rt,
    smiles_2D,
    inchikey_2D,
    score_biological,
    score_chemical,
    score_final,
    best_candidate_organism,
    consensus_1 = consensus_pat,
    consensus_2 = consensus_sup,
    consensus_3 = consensus_cla,
    consistency_1 = consistency_pat,
    consistency_2 = consistency_sup,
    consistency_3 = consistency_cla,
    best_candidate_1,
    best_candidate_2,
    best_candidate_3
  ) |>
  dplyr::mutate_all(list(~ y_as_na(x = ., y = ""))) |>
  dplyr::mutate(
    best_candidate_1 = if_else(
      condition = is.na(smiles_2D),
      true = "notAnnotated",
      false = best_candidate_1
    ),
    best_candidate_2 = if_else(
      condition = is.na(smiles_2D),
      true = paste(best_candidate_1, "notAnnotated"),
      false = best_candidate_2
    ),
    best_candidate_3 = if_else(
      condition = is.na(smiles_2D),
      true = paste(best_candidate_2, "notAnnotated"),
      false = best_candidate_3
    ),
    best_candidate_1 = if_else(
      condition = !is.na(smiles_2D) &
        is.na(best_candidate_1),
      true = "notClassified",
      false = best_candidate_1
    ),
    best_candidate_2 = if_else(
      condition = !is.na(smiles_2D) &
        is.na(best_candidate_2),
      true = paste(best_candidate_1, "notClassified"),
      false = best_candidate_2
    ),
    best_candidate_3 = if_else(
      condition = !is.na(smiles_2D) &
        is.na(best_candidate_3),
      true = paste(best_candidate_2, "notClassified"),
      false = best_candidate_3
    )
  )

## dirty TODO
colnames(feature_table) <-
  gsub(
    pattern = ".Peak.area",
    replacement = "",
    x = colnames(feature_table)
  )

log_debug(x = "removing \"row m/z\" and from \"row retention time\" columns")
feature_table <- feature_table |>
  dplyr::select(
    -"row m/z",
    -"row retention time",
    -"correlation group ID",
    -"annotation network number",
    -"best ion",
    -"auto MS2 verify",
    -"identified by n=",
    -"partners",
    -"neutral M mass",
    -"Unnamed: 64"
  ) |>
  tibble::column_to_rownames(var = "row ID")

top_n <- feature_table |>
  tibble::rownames_to_column() |>
  tidyr::gather(column, value, -rowname) |>
  dplyr::mutate(column = gsub(
    pattern = "^X",
    replacement = "",
    x = column
  )) |>
  dplyr::arrange(rowname, desc(value)) |>
  dplyr::filter(value >= INTENSITY)

top_m <- top_n |>
  dplyr::mutate(column = gsub(
    pattern = ".Peak.area",
    replacement = "",
    x = column
  )) |>
  dplyr::mutate(species = "Swertia chirayita") |>
  dplyr::select(
    feature_id = rowname,
    sample = column,
    intensity = value,
    species
  )

ms1_multiple <- ms1_best_candidate |>
  dplyr::left_join(top_m) |>
  dplyr::filter(!is.na(species)) |>
  dplyr::filter(intensity >= INTENSITY)

new_step <- ms1_multiple |>
  dplyr::mutate(
    rt_1 = as.numeric(rt),
    rt_2 = as.numeric(rt)
  ) |>
  data.table::data.table()

cat("setting joining keys \n")
data.table::setkey(df, rt_min, rt_max)
data.table::setkey(new_step, rt_1, rt_2)

cat("joining within given rt tolerance \n")
df_new_pre <- data.table::foverlaps(new_step, df)

df_new_with <- df_new_pre |>
  dplyr::rowwise() |>
  dplyr::filter(grepl(pattern = id, x = sample)) |>
  dplyr::mutate(
    mz_min = (1 - (1E-6 * PPM)) * as.numeric(mz),
    mz_max = (1 + (1E-6 * PPM)) * as.numeric(mz)
  ) |>
  dplyr::select(-rt_1, -rt_2) |>
  dplyr::filter(!is.na(peak_id))

df_new_without <- df_new_pre |>
  dplyr::filter(is.na(peak_id)) |>
  dplyr::filter(sample %in% df_new_with$sample)

df_peaks_samples <- list()

for (s in unique(df_new_with$id)) {
  df_new <- df_new_with |>
    dplyr::filter(id == s)

  df_peaks <- list()

  for (i in unique(df_new$peak_id)) {
    df_peak <- df_new |>
      dplyr::filter(peak_id == i) |>
      dplyr::arrange(desc(intensity))

    #' Silently sub-setting the object to speed-up analysis
    dda_data_min <-
      MSnbase::filterFile(
        dda_data,
        dda_data@phenoData@data$sampleNames[grepl(
          pattern = s,
          x = dda_data@phenoData@data$sampleNames
        )]
      ) |>
      MSnbase::filterRt(rt = c(
        min(df_peak$rt_min) * 60 - 10,
        max(df_peak$rt_max) * 60 + 10
      ))

    rtr <- df_peak |>
      mutate(rtmin = rt_min * 60, rtmax = rt_max * 60) |>
      select(rtmin, rtmax) |>
      as.matrix()

    mzr <- df_peak |>
      select(mzmin = mz_min, mzmax = mz_max) |>
      as.matrix()

    ms_chr <-
      chromatogram(
        object = dda_data_min,
        rt = rtr,
        mz = mzr
      )

    ## CAD part
    cad_time <- chromatograms_cad_improved[[1]][["time"]][which(
      chromatograms_cad_improved[[1]][["time"]] >= df_peak$rt_min[1] &
        chromatograms_cad_improved[[1]][["time"]] <= df_peak$rt_max[1]
    )] * 60
    cad_intensity <-
      chromatograms_cad_improved[[1]][["intensity"]][which(
        chromatograms_cad_improved[[1]][["time"]] >= df_peak$rt_min[1] &
          chromatograms_cad_improved[[1]][["time"]] <= df_peak$rt_max[1]
      )]

    cad_ready <-
      data.frame(intensity = cad_intensity, rtime = cad_time) |>
      dplyr::filter(!is.na(intensity)) |>
      dplyr::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
        min(intensity))) |>
      dplyr::filter(intensity >= 0.1) |>
      dplyr::mutate(rtime = (rtime - min(rtime)) / (max(rtime) -
        min(rtime))) |> # see https://github.com/sneumann/xcms/issues/593
      dplyr::arrange(rtime)

    cad_peak <-
      MSnbase::Chromatogram(intensity = cad_ready$intensity, rtime = cad_ready$rtime)

    ms_peaks <-
      parallel::mclapply(
        X = seq_along(ms_chr),
        FUN = extract_peak,
        mc.cores = 10
      )

    comparison_score <-
      parallel::mclapply(
        X = seq_along(ms_peaks),
        FUN = compare_peaks,
        mc.cores = 10
      )

    df_peak$comparison_score <- as.numeric(comparison_score)

    df_peak_new <- df_peak |>
      dplyr::filter(comparison_score >= PEAK_SIMILARITY_PREFILTER)

    df_peaks[[i]] <- df_peak_new
  }

  df_new_with_cor <- bind_rows(df_peaks)

  df_peaks_samples[[s]] <- df_new_with_cor
}

#' We got 2556 CAD peaks automatically detected
#' With peak similarity score > 0.6: 5313 features
df_new_with_cor_pre <- bind_rows(df_peaks_samples)

#' With peak similarity score > 0.7: 3421 features
df_new_with_cor_07 <- df_new_with_cor_pre |>
  dplyr::filter(comparison_score >= 0.7)

#' With peak similarity score > 0.8: 1928 features
df_new_with_cor_08 <- df_new_with_cor_pre |>
  dplyr::filter(comparison_score >= 0.8)

#' With peak similarity score > 0.75: 2598 features
df_new_with_cor_075 <- df_new_with_cor_pre |>
  dplyr::filter(comparison_score >= 0.75)

final_table_taxed <-
  prepare_hierarchy_preparation(dataframe = annotations) |>
  prepare_hierarchy()

final_table_taxed_with <-
  prepare_hierarchy(dataframe = df_new_with, detector = "ms")

final_table_taxed_without <-
  prepare_hierarchy(dataframe = df_new_without)

final_table_taxed_with_new <-
  prepare_hierarchy(
    dataframe = df_new_with,
    detector = "cad"
  )

final_table_taxed_with_new_cor_06 <-
  prepare_hierarchy(
    dataframe = df_new_with_cor_pre,
    detector = "cad"
  )
final_table_taxed_with_new_cor_07 <-
  prepare_hierarchy(
    dataframe = df_new_with_cor_07,
    detector = "cad"
  )
final_table_taxed_with_new_cor_08 <-
  prepare_hierarchy(
    dataframe = df_new_with_cor_08,
    detector = "cad"
  )
final_table_taxed_with_new_cor <-
  prepare_hierarchy(
    dataframe = df_new_with_cor_075,
    detector = "cad"
  )

samples <- prepare_plot(dataframe = final_table_taxed)
samples_with <- prepare_plot(dataframe = final_table_taxed_with)
samples_without <-
  prepare_plot(dataframe = final_table_taxed_without)
samples_with_new <-
  prepare_plot(dataframe = final_table_taxed_with_new)
samples_with_new_cor_06 <-
  prepare_plot(dataframe = final_table_taxed_with_new_cor_06)
samples_with_new_cor_07 <-
  prepare_plot(dataframe = final_table_taxed_with_new_cor_07)
samples_with_new_cor_08 <-
  prepare_plot(dataframe = final_table_taxed_with_new_cor_08)
samples_with_new_cor <-
  prepare_plot(dataframe = final_table_taxed_with_new_cor)

absolute <- plot_histograms(
  dataframe = samples,
  label = "Based on MS intensity only"
)

absolute_with <- plot_histograms(
  dataframe = samples_with,
  label = "MS intensity within CAD peak"
)

absolute_without <- plot_histograms(
  dataframe = samples_without,
  label = "MS intensity outside CAD peak"
)

absolute_with_new <- plot_histograms(
  dataframe = samples_with_new,
  label = "CAD intensity within CAD peak"
)

absolute_with_new_cor_06 <-
  plot_histograms(
    dataframe = samples_with_new_cor_06,
    label = "CAD intensity of corelated peaks within CAD peak"
  )
absolute_with_new_cor_07 <-
  plot_histograms(
    dataframe = samples_with_new_cor_07,
    label = "CAD intensity of corelated peaks within CAD peak"
  )
absolute_with_new_cor_08 <-
  plot_histograms(
    dataframe = samples_with_new_cor_08,
    label = "CAD intensity of corelated peaks within CAD peak"
  )
absolute_with_new_cor <-
  plot_histograms(
    dataframe = samples_with_new_cor,
    label = "CAD intensity of corelated peaks within CAD peak"
  )

combined <-
  absolute_with +
  absolute_without +
  absolute_with_new +
  absolute_with_new_cor

## specific sample exploration
plotly::plot_ly(
  data = final_table_taxed |>
    dplyr::filter(sample == "210619_AR_31_M_36_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "sunburst",
  branchvalues = "total"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_with |>
    dplyr::filter(sample == "210619_AR_31_M_36_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "sunburst",
  branchvalues = "total"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_without |>
    dplyr::filter(sample == "210619_AR_31_M_36_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "sunburst",
  branchvalues = "total"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_with_new |>
    dplyr::filter(sample == "210619_AR_31_M_36_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "sunburst",
  branchvalues = "total"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_with_new_cor |>
    dplyr::filter(sample == "210619_AR_31_M_36_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "sunburst",
  branchvalues = "total"
) |>
  plotly::layout(colorway = sunburst_colors)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
