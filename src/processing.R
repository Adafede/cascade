start <- Sys.time()

library(package = baseline, quietly = TRUE)
library(package = data.table, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = docopt, quietly = TRUE)
library(package = microshades, quietly = TRUE)
library(package = MSnbase, quietly = TRUE)
library(package = nucleR, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = pracma, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = readr, quietly = TRUE)

source(file = "R/colors.R")
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

FOURRIER_COMPONENTS <- 0.05 ## of total

TIME_MIN <- 0 ## Min

TIME_MAX <- 32 ## Min

FREQUENCY <- 20 ## Hz

RESAMPLE <- 1 / 10 ## points

CAD_SHIFT <- 0.055 ## minutes

PDA_SHIFT <- 0.090 ## minutes

ESTIMATED_SOLUBLITIY_LIMIT <- 58

WORKERS <- 10

# adapted from Excel sheet from paper shortDOI: 10/ghmvhz
sigma <- 0.05
k2 <- sigma / 250 # 30
k4 <- sigma / 1250000 # 200
smoothing_width <- 5
baseline_adjust <- 0

toyset <- "~/data/20210701_10043/test"

future::plan(strategy = future::multiprocess(workers = WORKERS))

files <- list.files(
  path = toyset,
  pattern = ".mzML.gz",
  full.names = TRUE,
  recursive = TRUE
)

names <- list.files(
  path = toyset,
  pattern = ".mzML.gz",
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML.gz",
    replacement = "",
    fixed = TRUE
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

chromatograms_bpi <- chromatograms_all[c(TRUE, FALSE, FALSE)]

chromatograms_pda <- chromatograms_all[c(FALSE, TRUE, FALSE)]

chromatograms_cad <- chromatograms_all[c(FALSE, FALSE, TRUE)]

chromatograms_bpi_improved <- list()

for (i in seq_along(1:length(chromatograms_bpi))) {
  chromatograms_bpi_improved[[i]] <-
    improve_signal(df = chromatograms_bpi[[i]] |>
      dplyr::select(time, intensity = BasePeak_0))
}

names(chromatograms_bpi_improved) <- names

bpis_improved <-
  dplyr::bind_rows(chromatograms_bpi_improved, .id = "id") |>
  dplyr::mutate(time = time)

chromatograms_pda_improved <- list()

for (i in seq_along(1:length(chromatograms_pda))) {
  chromatograms_pda_improved[[i]] <-
    improve_signal(df = chromatograms_pda[[i]] |>
      dplyr::select(time, intensity = PDA.1_TotalAbsorbance_0))
}

names(chromatograms_pda_improved) <- names

pdas_improved <-
  dplyr::bind_rows(chromatograms_pda_improved, .id = "id") |>
  dplyr::mutate(time = time + PDA_SHIFT)

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

bpi_plot <- plotly::plot_ly(
  data = bpis_improved,
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
      text = "MS (POS)",
      showarrow = FALSE
    )
  )

pda_plot <- plotly::plot_ly(
  data = pdas_improved,
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
      text = "UV (200-500nm)",
      showarrow = FALSE
    )
  )

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

comparison <-
  plotly::subplot(bpi_plot,
    pda_plot,
    cad_plot,
    nrows = 3,
    shareX = TRUE
  )

comparison

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

comparison_2 <-
  plotly::subplot(bpi_plot,
    pda_plot,
    new_new,
    nrows = 3,
    shareX = TRUE
  )

comparison_2

# Not run:
plot(cads_baselined$intensity, type = "l", col = "navy")
grid()
x <-
  findpeaks(
    cads_baselined$intensity,
    npeaks = 2000,
    threshold = 0.01,
    sortstr = TRUE
  )
points(x[, 2], x[, 1], pch = 20, col = "maroon") ## End(Not run)

peaks <- data.frame(x) |>
  mutate(
    peak_id = row_number(),
    peak_max = X1,
    rt_apex = cads_baselined$time[X2],
    rt_min = cads_baselined$time[X3],
    rt_max = cads_baselined$time[X4]
  ) |>
  data.table()

cads_baselined <- cads_baselined |>
  mutate(rt_1 = time, rt_2 = time) |>
  data.table()

cat("setting joining keys \n")
setkey(peaks, rt_min, rt_max)
setkey(cads_baselined, rt_1, rt_2)

cat("joining within given rt tolerance \n")
df <- foverlaps(peaks, cads_baselined) |>
  group_by(peak_id) |>
  mutate(integral = sum(intensity)) |>
  ungroup() |>
  select(peak_id, peak_max, rt_apex, rt_min, rt_max, integral) |>
  distinct() |>
  mutate(integral = integral / sum(integral)) |>
  filter(integral >= 0.01) |>
  data.table()

annotations <-
  readr::read_delim(file = "~/git/tima-r/data/processed/210718_163140/yourFinalFile.tsv.gz")

clean_xanthones <- TRUE

ms1_best_candidate <- annotations |>
  dplyr::mutate_all(list(~ gsub(
    pattern = "\\|.*",
    replacement = "",
    x = .x
  ))) |>
  splitstackshape::cSplit("best_candidate", sep = "§") |>
  dplyr::distinct(
    feature_id,
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
feature_table <-
  read_features(id = "97d7c50031a84b9ba2747e883a5907cd")

colnames(feature_table) <-
  gsub(
    pattern = ".Peak.area",
    replacement = "",
    x = colnames(feature_table)
  )

log_debug(x = "removing \"row m/z\" and from \"row retention time\" columns")
feature_table <- feature_table %>%
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
  dplyr::arrange(rowname, desc(value))

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
  ## add this step
  dplyr::filter(!is.na(species)) |>
  dplyr::filter(intensity != 0)

new_step <- ms1_multiple |>
  mutate(
    rt_1 = as.numeric(rt),
    rt_2 = as.numeric(rt)
  ) |>
  data.table()

cat("setting joining keys \n")
setkey(df, rt_min, rt_max)
setkey(new_step, rt_1, rt_2)

cat("joining within given rt tolerance \n")
df_new <- foverlaps(new_step, df)

df_new_with <- df_new |>
  filter(!is.na(peak_id))

df_new_without <- df_new |>
  filter(is.na(peak_id))

final_table_taxed <- annotations |>
  prepare_hierarchy_preparation() |>
  prepare_hierarchy() |>
  dplyr::mutate(species = "Swertia chirayita")

final_table_taxed_with <- df_new_with |>
  prepare_hierarchy() |>
  dplyr::mutate(species = "Swertia chirayita")

final_table_taxed_without <- df_new_without |>
  prepare_hierarchy() |>
  dplyr::mutate(species = "Swertia chirayita")

nice_colors <- rev(
  list(
    microshades_palette("micro_cvd_green", lightest = FALSE),
    microshades_palette("micro_cvd_orange", lightest = FALSE),
    microshades_palette("micro_cvd_blue", lightest = FALSE),
    microshades_palette("micro_cvd_turquoise", lightest = FALSE),
    microshades_palette("micro_cvd_purple", lightest = FALSE),
    microshades_palette("micro_cvd_gray", lightest = FALSE),
    microshades_palette("micro_orange", lightest = FALSE),
    microshades_palette("micro_purple", lightest = FALSE)
  )
)

sunburst_colors <- character()

for (i in seq_len(length(nice_colors))) {
  sunburst_colors[i] <- rev(nice_colors)[[i]][5]
}

samples <- prepare_plot(dataframe = final_table_taxed)
samples_with <- prepare_plot(dataframe = final_table_taxed_with)
samples_without <-
  prepare_plot(dataframe = final_table_taxed_without)

absolute <-
  plot_histograms(dataframe = samples, label = "Based on MS intensity only")

absolute_with <-
  plot_histograms(dataframe = samples_with, label = "MS intensity within CAD peak")

absolute_without <-
  plot_histograms(dataframe = samples_without, label = "MS intensity outside CAD peak")

ggpubr::ggarrange(
  absolute,
  absolute_with,
  absolute_without,
  nrow = 3,
  align = "v",
  common.legend = TRUE,
  legend = "right"
)

## specific sample exploration
plotly::plot_ly(
  data = final_table_taxed |>
    dplyr::filter(species == "Swertia chirayita") |>
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
    dplyr::filter(species == "Swertia chirayita") |>
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
    dplyr::filter(species == "Swertia chirayita") |>
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
