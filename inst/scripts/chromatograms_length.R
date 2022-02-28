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
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
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
WORKERS <- 10

#' Parameters for LC alignment
TIME_MIN <- 0 ## minutes
TIME_MAX <- 999 ## minutes
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
RT_TOL <- 0.1
PPM <- 10

#' Parameters for annotation
CONFIDENCE_SCORE_MIN <- 0.5

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

chromatograms_cad_improved <-
  improve_signals_progress(chromatograms_cad)

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
                               nrows = 3) |>
  plotly::layout(yaxis = list(range = c(0, 0.2)),
                 yaxis2 = list(range = c(0, 0.2)),
                 yaxis3 = list(range = c(0, 0.2)))
new_new_new

peaks_cad <- peaks_progress(chromatograms_cad_baselined)

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
  dplyr::arrange(rowname, dplyr::desc(value)) |>
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
  dplyr::ungroup() |>
  dplyr::select(-rt_1, -rt_2) |>
  dplyr::filter(!is.na(peak_id)) |>
  make_confident(score = CONFIDENCE_SCORE_MIN)

df_new_without <- df_new_pre |>
  dplyr::filter(is.na(peak_id)) |>
  dplyr::filter(sample %in% df_new_with$sample) |>
  make_confident(score = CONFIDENCE_SCORE_MIN)

df_peaks_samples <- list()

for (s in unique(df_new_with$id)) {
  df_new <- df_new_with |>
    dplyr::filter(id == s)

  df_peaks <- list()

  for (i in unique(df_new$peak_id)) {
    df_peak <- df_new |>
      dplyr::filter(peak_id == i) |>
      dplyr::arrange(dplyr::desc(intensity))

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
        (min(df_peak$rt_min) + CAD_SHIFT - RT_TOL) * 60,
        (max(df_peak$rt_max) + CAD_SHIFT + RT_TOL) * 60
      )) |>
      MSnbase::filterMz(mz = c(
        min(df_peak$mz_min),
        max(df_peak$mz_max)
      ))

    ## CAD part
    cad_ready <- chromatograms_cad_improved[[s]] |>
      dplyr::filter(time >= df_peak$rt_min[1] &
        time <= df_peak$rt_max[1]) |>
      dplyr::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
        min(intensity))) |>
      dplyr::filter(intensity >= 0.1)

    cad_peak <-
      MSnbase::Chromatogram(
        intensity = cad_ready$intensity,
        rtime = (cad_ready$time - min(cad_ready$time)) /
          (max(cad_ready$time) - min(cad_ready$time))
      ) #' see https://github.com/sneumann/xcms/issues/593

    ms_peaks <-
      parallel::mclapply(
        X = seq_along(1:nrow(df_peak)),
        FUN = extract_peak
      )

    comparison_score <-
      if (length(ms_peaks) != 0) {
        parallel::mclapply(
          X = seq_along(ms_peaks),
          FUN = compare_peaks
        )
      } else {
        rep(list(0), nrow(df_peak))
      }

    df_peak$comparison_score <- as.numeric(comparison_score)

    df_peak_new <- df_peak

    df_peaks[[i]] <- df_peak_new
  }

  df_new_with_cor <- dplyr::bind_rows(df_peaks)

  df_peaks_samples[[s]] <- df_new_with_cor
}

df_peaks_samples_full <- dplyr::bind_rows(df_peaks_samples)

#' We got 2556 CAD peaks automatically detected
#' With peak similarity score > 0.6: 5313 features
df_new_with_cor_pre <- df_peaks_samples_full |>
  dplyr::filter(comparison_score >= PEAK_SIMILARITY_PREFILTER) #' TODO check negative values

df_new_with_cor_pre_taxo <- df_new_with_cor_pre |>
  dplyr::rowwise() |>
  dplyr::mutate(taxo = ifelse(
    test = best_candidate_organism %in% species,
    yes = 1,
    no = 0
  )) |>
  dplyr::group_by(id, peak_id) |>
  dplyr::mutate(sum = sum(taxo)) |>
  filter(taxo == 1 | sum == 0)

#' dirty annotation change after computation of peak comparison
# df_new_with_cor_pre <- df_new_with_cor_pre |>
#   dplyr::select(
#     id,
#     peak_id,
#     peak_max,
#     rt_apex,
#     rt_min,
#     rt_max,
#     integral,
#     feature_id,
#     sample,
#     intensity,
#     species,
#     mz_min,
#     mz_max,
#     comparison_score
#   ) |>
#   dplyr::left_join(ms1_multiple |>
#                      mutate(feature_id = as.character(feature_id))) |>
#   make_confident(score = CONFIDENCE_SCORE_MIN)

#' With peak similarity score > 0.7: 3421 features
df_new_with_cor_07 <- df_new_with_cor_pre_taxo |>
  dplyr::filter(comparison_score >= 0.7)

#' With peak similarity score > 0.8: 1928 features
df_new_with_cor_08 <- df_new_with_cor_pre_taxo |>
  dplyr::filter(comparison_score >= 0.8)

#' With peak similarity score > 0.75: 2598 features
df_new_with_cor <- df_new_with_cor_pre_taxo |>
  dplyr::filter(comparison_score >= 0.75)

#' TODO fix it
# final_table_taxed <-
#   prepare_hierarchy_preparation(dataframe = annotations) |>
#   prepare_hierarchy()

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
    dataframe = df_new_with_cor,
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

bound <- dplyr::bind_rows(
  df_new_with |> dplyr::mutate(species = "absolute_with"),
  df_new_without |> dplyr::mutate(species = "absolute_without"),
  df_new_with |> dplyr::mutate(intensity = integral, species = "absolute_with_new"),
  df_new_with_cor_075 |> dplyr::mutate(intensity = integral, species = "absolute_with_new_cor")
)

bound_ready <- bound |>
  prepare_hierarchy(rescale = TRUE) |>
  prepare_plot()

absolute <- plot_histograms(
  dataframe = samples,
  label = "Based on MS intensity only",
  xlab = FALSE
)

absolute_with <- plot_histograms(
  dataframe = samples_with,
  # dataframe = bound_ready[bound_ready$species == "absolute_with", ] |>
  #   dplyr::mutate_at(c("ids", "sample", "color"),as.character) |>
  #   dplyr::arrange(sample,parents,ids) |>
  #   dplyr::mutate_at(c("ids", "sample", "color"),as.factor),
  label = "MS intensity within CAD peak",
  xlab = FALSE
)

absolute_without <- plot_histograms(
  dataframe = samples_without,
  # dataframe = bound_ready[bound_ready$species == "absolute_without", ] |>
  #   dplyr::mutate_at(c("ids", "sample", "color"),as.character) |>
  #   dplyr::arrange(sample,parents,ids) |>
  #   dplyr::mutate_at(c("ids", "sample", "color"),as.factor),
  label = "MS intensity outside CAD peak",
  xlab = FALSE
)

absolute_with_new <- plot_histograms(
  dataframe = samples_with_new,
  # dataframe = bound_ready[bound_ready$species == "absolute_with_new", ] |>
  #   dplyr::mutate_at(c("ids", "sample", "color"),as.character) |>
  #   dplyr::arrange(sample,parents,ids) |>
  #   dplyr::mutate_at(c("ids", "sample", "color"),as.factor),
  label = "CAD intensity within CAD peak",
  xlab = FALSE
)

absolute_with_new_cor_06 <-
  plot_histograms(
    dataframe = samples_with_new_cor_06,
    label = "CAD intensity of corelated peaks within CAD peak",
    xlab = FALSE
  )
absolute_with_new_cor_07 <-
  plot_histograms(
    dataframe = samples_with_new_cor_07,
    label = "CAD intensity of corelated peaks within CAD peak",
    xlab = FALSE
  )
absolute_with_new_cor_08 <-
  plot_histograms(
    dataframe = samples_with_new_cor_08,
    label = "CAD intensity of corelated peaks within CAD peak",
    xlab = FALSE
  )
absolute_with_new_cor <-
  plot_histograms(
    dataframe = samples_with_new_cor,
    # dataframe = bound_ready[bound_ready$species == "absolute_with_new_cor", ] |>
    #   dplyr::mutate_at(c("ids", "sample", "color"),as.character) |>
    #   dplyr::arrange(sample,parents,ids) |>
    #   dplyr::mutate_at(c("ids", "sample", "color"),as.factor),
    label = "CAD intensity of corelated peaks within CAD peak",
    xlab = FALSE
  )

combined <-
  absolute_with +
  absolute_without +
  absolute_with_new +
  absolute_with_new_cor

combined_2 <-
  absolute_with_new_cor_06 +
  absolute_with_new_cor_07 +
  absolute_with_new_cor_08 +
  absolute_with_new_cor

## specific sample exploration
plotly::plot_ly(
  data = final_table_taxed |>
    dplyr::filter(sample == "210619_AR_25_M_30_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "treemap",
  branchvalues = "total",
  textinfo = "label+value+percent parent+percent root"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_with |>
    dplyr::filter(sample == "210619_AR_25_M_30_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "treemap",
  branchvalues = "total",
  textinfo = "label+value+percent parent+percent root"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_without |>
    dplyr::filter(sample == "210619_AR_25_M_30_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "treemap",
  branchvalues = "total",
  textinfo = "label+value+percent parent+percent root"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_with_new |>
    dplyr::filter(sample == "210619_AR_25_M_30_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "treemap",
  branchvalues = "total",
  textinfo = "label+value+percent parent+percent root"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_with_new_cor_07 |>
    dplyr::filter(sample == "210619_AR_25_M_30_01"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "treemap",
  branchvalues = "total",
  textinfo = "label+value+percent parent+percent root"
) |>
  plotly::layout(colorway = sunburst_colors)

# plotly::plot_ly(
#   table_new  |> #' from internal prepare_hierarchy()
#     dplyr::filter(!grepl(pattern = " Other-",
#                          x = parents,
#                          fixed = TRUE)) |>
#     dplyr::mutate_at(c("ids", "sample"), as.character) |>
#     dplyr::arrange(sample, parents, ids) |>
#     dplyr::mutate_at(c("ids", "sample"), as.factor) |>
#     dplyr::left_join(df_new_with_cor_075 |> distinct(feature_id, rt_apex)),
#   x = ~ rt_apex,
#   y = ~ intensity,
#   type = "bar",
#   color = ~ parents
# ) |>
#   plotly::layout(barmode = "stack")

#' Work in progress
#' Add some metadata per peak
df_meta <- df_new_with_cor_pre_taxo |>
  dplyr::arrange(desc(intensity)) |>
  dplyr::distinct(
    id,
    peak_id,
    integral,
    feature_id,
    rt,
    mz,
    smiles_2D,
    inchikey_2D,
    best_candidate_1,
    best_candidate_2,
    best_candidate_3,
    score_biological,
    score_chemical,
    score_final,
    #' add consensus
    sample,
    species
  ) |>
  dplyr::group_by(id, peak_id) |>
  dplyr::distinct(feature_id,
    # smiles_2D,
    # inchikey_2D,
    # best_candidate_1,
    # best_candidate_2,
    # best_candidate_3,
    .keep_all = TRUE
  ) |>
  dplyr::add_count(name = "featuresPerPeak") |>
  dplyr::distinct(smiles_2D,
    inchikey_2D,
    # best_candidate_1,
    # best_candidate_2,
    # best_candidate_3,
    .keep_all = TRUE
  ) |>
  dplyr::add_count(name = "structuresPerPeak") |>
  dplyr::distinct(best_candidate_1,
    best_candidate_2,
    best_candidate_3,
    .keep_all = TRUE
  ) |>
  dplyr::add_count(name = "chemicalClassesPerPeak") |>
  dplyr::distinct(best_candidate_1,
    best_candidate_2,
    .keep_all = TRUE
  ) |>
  dplyr::add_count(name = "chemicalSuperclassesPerPeak") |>
  dplyr::distinct(best_candidate_1,
    .keep_all = TRUE
  ) |>
  dplyr::add_count(name = "chemicalPathwaysPerPeak") |>
  dplyr::ungroup() |>
  dplyr::distinct(
    id,
    peak_id,
    featuresPerPeak,
    structuresPerPeak,
    chemicalClassesPerPeak,
    chemicalSuperclassesPerPeak,
    chemicalPathwaysPerPeak
  )

test_features <- df_meta |>
  dplyr::group_by(featuresPerPeak) |>
  dplyr::count()
test_structures <- df_meta |>
  dplyr::group_by(structuresPerPeak) |>
  dplyr::count()
test_classes <- df_meta |>
  dplyr::group_by(chemicalClassesPerPeak) |>
  dplyr::count()
test_superclasses <- df_meta |>
  dplyr::group_by(chemicalSuperclassesPerPeak) |>
  dplyr::count()
test_pathways <- df_meta |>
  dplyr::group_by(chemicalPathwaysPerPeak) |>
  dplyr::count()

plotly::plot_ly() |>
  plotly::add_pie(
    data = test_features,
    name = "Features",
    labels = ~featuresPerPeak,
    values = ~n,
    sort = FALSE,
    type = "pie",
    textposition = "inside",
    domain = list(row = 0, column = 0)
  ) |>
  plotly::add_pie(
    data = test_structures,
    name = "Structures",
    labels = ~structuresPerPeak,
    values = ~n,
    type = "pie",
    textposition = "inside",
    domain = list(row = 0, column = 1)
  ) |>
  plotly::add_pie(
    data = test_classes,
    name = "Classes",
    labels = ~chemicalClassesPerPeak,
    values = ~n,
    type = "pie",
    textposition = "inside",
    domain = list(row = 0, column = 2)
  ) |>
  plotly::add_pie(
    data = test_superclasses,
    name = "Superclasses",
    labels = ~chemicalSuperclassesPerPeak,
    values = ~n,
    type = "pie",
    textposition = "inside",
    domain = list(row = 1, column = 0)
  ) |>
  plotly::add_pie(
    data = test_pathways,
    name = "Pathways",
    labels = ~chemicalPathwaysPerPeak,
    values = ~n,
    type = "pie",
    textposition = "inside",
    domain = list(row = 1, column = 1)
  ) |>
  plotly::layout(
    title = "Peak analysis \n Features > Structures > Chemical classes > \n Superclasses > Pathways",
    colorway = viridis::cividis(max(
      nrow(test_features),
      nrow(test_structures),
      nrow(test_classes),
      nrow(test_superclasses),
      nrow(test_pathways)
    )),
    grid = list(rows = 2, columns = 3)
  )

end <- Sys.time()

log_debug("Script finished in", format(end - start))
