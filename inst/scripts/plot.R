start <- Sys.time()

library(package = data.table, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = parallel, quietly = TRUE)
library(package = patchwork, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = readr, quietly = TRUE)

source(file = "R/check_export_dir.R")
source(file = "R/colors.R")
source(file = "R/get_gnps.R")
source(file = "R/get_params.R")
source(file = "R/log_debug.R")
source(file = "R/make_confident.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/plot_histograms.R")
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

#' Paths
ANNOTATIONS <-
  "~/git/tima-r/inst/extdata/processed/220208_172733/20220208_10043.tsv.gz"
EXPORT_DIR <- "~/git/cascade/data/interim"
EXPORT_FILE <- "peaks_informed.tsv.gz"

#' Parameters related to MS/CAD
INTENSITY_MS_MIN <- 1E5
PEAK_SIMILARITY <- 0.9
PEAK_SIMILARITY_PREFILTER <- 0.6
RT_TOL <- 0.1
PPM <- 10
AREA_MIN <- 0.005

#' Parameters for annotation
CONFIDENCE_SCORE_MIN <- 0.5

log_debug(x = "loading compared peaks")
compared_peaks <-
  readr::read_delim(file = file.path(EXPORT_DIR, EXPORT_FILE))

log_debug(x = "loading annotations")
annotations <- readr::read_delim(file = ANNOTATIONS)

log_debug(x = "keeping best annotations only")
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

log_debug(x = "adding metadata dirtily for now")
candidates_metadata <- ms1_best_candidate |>
  dplyr::mutate(species = "Swertia chirayita") |>
  dplyr::mutate(feature_id = as.numeric(feature_id))

log_debug(x = "keeping only candidates above desired threshold")
candidates_confident <- candidates_metadata |>
  make_confident(score = CONFIDENCE_SCORE_MIN)

log_debug(x = "joining compared peaks and candidates")
df_peaks_samples_full <- compared_peaks |>
  dplyr::left_join(candidates_confident)

log_debug(x = "temporary fix") #' TODO
df_peaks_samples_full <- df_peaks_samples_full |>
  dplyr::mutate(
    id = sample,
    integral = peak_area,
    intensity = feature_area
  )

log_debug(x = "keeping peaks similarities above desired (pre-)threshold only")
df_new_with_cor_pre <- df_peaks_samples_full |>
  dplyr::filter(comparison_score >= PEAK_SIMILARITY_PREFILTER) #' TODO check negative values

log_debug(x = "keeping multiple features only if none was reported in the sepcies")
df_new_with_cor_pre_taxo <- df_new_with_cor_pre |>
  dplyr::rowwise() |>
  dplyr::mutate(taxo = ifelse(
    test = best_candidate_organism %in% species,
    yes = 1,
    no = 0
  )) |>
  dplyr::group_by(sample, peak_id) |> #' TODO switch to ID if needed
  dplyr::mutate(sum = sum(taxo)) |>
  filter(taxo == 1 | sum == 0)

log_debug(x = "keeping peaks similarities with score above ...")
log_debug(x = "... 0.7")
df_new_with_cor_07 <- df_new_with_cor_pre_taxo |>
  dplyr::filter(comparison_score >= 0.7)

log_debug(x = "... 0.8")
df_new_with_cor_08 <- df_new_with_cor_pre_taxo |>
  dplyr::filter(comparison_score >= 0.8)

log_debug(x = "... 0.75")
df_new_with_cor <- df_new_with_cor_pre_taxo |>
  dplyr::filter(comparison_score >= 0.75)

#' TODO fix it
# final_table_taxed <-
#   prepare_hierarchy_preparation(dataframe = annotations) |>
#   prepare_hierarchy()
# final_table_taxed_with <-
#   prepare_hierarchy(dataframe = df_new_with, detector = "ms")
#
# final_table_taxed_without <-
#   prepare_hierarchy(dataframe = df_new_without)
#
# final_table_taxed_with_new <-
#   prepare_hierarchy(
#     dataframe = df_new_with,
#     detector = "cad"
#   )

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

#' TODO fix it
# samples <- prepare_plot(dataframe = final_table_taxed)
# samples_with <- prepare_plot(dataframe = final_table_taxed_with)
# samples_without <-
#   prepare_plot(dataframe = final_table_taxed_without)
# samples_with_new <-
#   prepare_plot(dataframe = final_table_taxed_with_new)
samples_with_new_cor_06 <-
  prepare_plot(dataframe = final_table_taxed_with_new_cor_06)
samples_with_new_cor_07 <-
  prepare_plot(dataframe = final_table_taxed_with_new_cor_07)
samples_with_new_cor_08 <-
  prepare_plot(dataframe = final_table_taxed_with_new_cor_08)
samples_with_new_cor <-
  prepare_plot(dataframe = final_table_taxed_with_new_cor)

# bound <- dplyr::bind_rows(
#   df_new_with |> dplyr::mutate(species = "absolute_with"),
#   df_new_without |> dplyr::mutate(species = "absolute_without"),
#   df_new_with |> dplyr::mutate(intensity = integral, species = "absolute_with_new"),
#   df_new_with_cor |> dplyr::mutate(intensity = integral, species = "absolute_with_new_cor")
# )

# bound_ready <- bound |>
#   prepare_hierarchy(rescale = TRUE) |>
#   prepare_plot()

# absolute <- plot_histograms(
#   dataframe = samples,
#   label = "Based on MS intensity only",
#   xlab = FALSE
# )

# absolute_with <- plot_histograms(
#   dataframe = samples_with,
#   # dataframe = bound_ready[bound_ready$species == "absolute_with", ] |>
#   #   dplyr::mutate_at(c("ids", "sample", "color"),as.character) |>
#   #   dplyr::arrange(sample,parents,ids) |>
#   #   dplyr::mutate_at(c("ids", "sample", "color"),as.factor),
#   label = "MS intensity within CAD peak",
#   xlab = FALSE
# )

# absolute_without <- plot_histograms(
#   dataframe = samples_without,
#   # dataframe = bound_ready[bound_ready$species == "absolute_without", ] |>
#   #   dplyr::mutate_at(c("ids", "sample", "color"),as.character) |>
#   #   dplyr::arrange(sample,parents,ids) |>
#   #   dplyr::mutate_at(c("ids", "sample", "color"),as.factor),
#   label = "MS intensity outside CAD peak",
#   xlab = FALSE
# )

# absolute_with_new <- plot_histograms(
#   dataframe = samples_with_new,
#   # dataframe = bound_ready[bound_ready$species == "absolute_with_new", ] |>
#   #   dplyr::mutate_at(c("ids", "sample", "color"),as.character) |>
#   #   dplyr::arrange(sample,parents,ids) |>
#   #   dplyr::mutate_at(c("ids", "sample", "color"),as.factor),
#   label = "CAD intensity within CAD peak",
#   xlab = FALSE
# )

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

# combined <-
#   absolute_with +
#   absolute_without +
#   absolute_with_new +
#   absolute_with_new_cor

combined_2 <-
  absolute_with_new_cor_06 +
  absolute_with_new_cor_07 +
  absolute_with_new_cor_08 +
  absolute_with_new_cor

## specific sample exploration
# plotly::plot_ly(
#   data = final_table_taxed |>
#     dplyr::filter(sample == "210619_AR_25_M_30_01"),
#   ids = ~ids,
#   labels = ~labels,
#   parents = ~parents,
#   values = ~values,
#   maxdepth = 3,
#   type = "treemap",
#   branchvalues = "total",
#   textinfo = "label+value+percent parent+percent root"
# ) |>
#   plotly::layout(colorway = sunburst_colors)

# plotly::plot_ly(
#   data = final_table_taxed_with |>
#     dplyr::filter(sample == "210619_AR_25_M_30_01"),
#   ids = ~ids,
#   labels = ~labels,
#   parents = ~parents,
#   values = ~values,
#   maxdepth = 3,
#   type = "treemap",
#   branchvalues = "total",
#   textinfo = "label+value+percent parent+percent root"
# ) |>
#   plotly::layout(colorway = sunburst_colors)

# plotly::plot_ly(
#   data = final_table_taxed_without |>
#     dplyr::filter(sample == "210619_AR_25_M_30_01"),
#   ids = ~ids,
#   labels = ~labels,
#   parents = ~parents,
#   values = ~values,
#   maxdepth = 3,
#   type = "treemap",
#   branchvalues = "total",
#   textinfo = "label+value+percent parent+percent root"
# ) |>
#   plotly::layout(colorway = sunburst_colors)

# plotly::plot_ly(
#   data = final_table_taxed_with_new |>
#     dplyr::filter(sample == "210619_AR_25_M_30_01"),
#   ids = ~ids,
#   labels = ~labels,
#   parents = ~parents,
#   values = ~values,
#   maxdepth = 3,
#   type = "treemap",
#   branchvalues = "total",
#   textinfo = "label+value+percent parent+percent root"
# ) |>
#   plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = final_table_taxed_with_new_cor_07,
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

plotly::plot_ly(
  data = df_new_with_cor,
  x = ~peak_rt_apex,
  y = ~peak_area,
  color = ~best_candidate_1,
  colors = rev(paired),
  type = "bar"
) |>
  plotly::layout(
    barmode = "stack"
  )

end <- Sys.time()

log_debug("Script finished in", format(end - start))
