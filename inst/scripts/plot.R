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
source(file = "R/make_other.R")
source(file = "R/no_other.R")
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
  file.path(
    paths$inst$extdata$interim$annotations$path,
    params$annotation$tool,
    params$filename$annotation$tima
  )

EXPORT_DIR <- paths$inst$extdata$interim$peaks
EXPORT_FILE_CAD <- list.files(
  path = EXPORT_DIR,
  pattern = paste(params$filename$mzml, "featuresInformed_cad", sep = "_"),
  full.names = FALSE,
  recursive = FALSE
)
EXPORT_FILE_CAD_2 <- list.files(
  path = EXPORT_DIR,
  pattern = paste(params$filename$mzml, "featuresNotInformed_cad", sep = "_"),
  full.names = FALSE,
  recursive = FALSE
)
EXPORT_FILE_PDA <- list.files(
  path = EXPORT_DIR,
  pattern = paste(params$filename$mzml, "featuresInformed_pda", sep = "_"),
  full.names = FALSE,
  recursive = FALSE
)
EXPORT_FILE_PDA_2 <- list.files(
  path = EXPORT_DIR,
  pattern = paste(params$filename$mzml, "featuresNotInformed_pda", sep = "_"),
  full.names = FALSE,
  recursive = FALSE
)

#' Parameters related to MS/CAD
INTENSITY_MS_MIN <- params$chromato$intensity$ms1$min
PEAK_SIMILARITY <- params$chromato$peak$similarity$filter
PEAK_SIMILARITY_PREFILTER <-
  params$chromato$peak$similarity$prefilter
RT_TOL <- params$chromato$peak$tolerance$rt
PPM <- params$chromato$peak$tolerance$ppm
AREA_MIN <- params$chromato$peak$area$min

#' Parameters for annotation
CONFIDENCE_SCORE_MIN <- params$annotation$confidence$min

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

if (params$signal$detector$cad == TRUE) {
  log_debug(x = "loading compared peaks")
  compared_peaks_cad <-
    readr::read_delim(file = file.path(EXPORT_DIR, EXPORT_FILE_CAD))

  outside_peaks_cad <-
    readr::read_delim(file = file.path(EXPORT_DIR, EXPORT_FILE_CAD_2))

  log_debug(x = "joining compared peaks and candidates")
  df_peaks_samples_full_cad <- compared_peaks_cad |>
    dplyr::left_join(candidates_confident)

  df_peaks_samples_min_cad <- compared_peaks_cad |>
    dplyr::left_join(candidates_confident)

  log_debug(x = "temporary fix") #' TODO
  df_peaks_samples_full_cad <- df_peaks_samples_full_cad |>
    dplyr::mutate(
      id = sample,
      integral = peak_area,
      intensity = feature_area
    )
  df_peaks_samples_min_cad <- df_peaks_samples_min_cad |>
    dplyr::mutate(
      id = sample,
      integral = peak_area,
      intensity = feature_area,
      peak_rt_apex = rt,
      peak_area = 1
    )

  log_debug(x = "keeping peaks similarities above desired (pre-)threshold only")
  df_new_with_cor_pre_cad <- df_peaks_samples_full_cad |>
    dplyr::filter(comparison_score >= PEAK_SIMILARITY_PREFILTER) #' TODO check negative values

  log_debug(x = "keeping multiple features only if none was reported in the sepcies")
  df_new_with_cor_pre_taxo_cad <- df_new_with_cor_pre_cad |>
    dplyr::rowwise() |>
    dplyr::mutate(taxo = ifelse(
      test = grepl(
        pattern = best_candidate_organism,
        x = species
      ),
      yes = 1,
      no = 0
    )) |>
    dplyr::mutate(taxo = ifelse(
      test = is.na(taxo),
      yes = 0,
      no = taxo
    )) |>
    dplyr::mutate(taxo_2 = ifelse(
      test = best_candidate_organism %in% species,
      yes = 1,
      no = 0
    )) |>
    dplyr::group_by(sample, peak_id) |> #' TODO switch to ID if needed
    dplyr::mutate(sum = sum(taxo)) |>
    dplyr::mutate(sum_2 = sum(taxo_2)) |>
    dplyr::filter(taxo_2 == 1 |
      sum_2 == 0 |
      sum == 0) |> #' TODO decide if genus
    dplyr::ungroup()

  log_debug(x = "keeping peaks similarities with score above", PEAK_SIMILARITY)
  df_new_with_cor_cad <- df_new_with_cor_pre_taxo_cad |>
    dplyr::filter(comparison_score >= PEAK_SIMILARITY)

  log_debug(x = "plotting histograms")
  #' TODO harmonize 'others' among minor and major
  df_histogram_ready_cad <- df_new_with_cor_pre_taxo_cad |>
    make_other() |>
    prepare_plot_2()
  df_histogram_outside_ready_cad <- df_peaks_samples_min_cad |>
    make_other() |>
    prepare_plot_2()

  df_histogram_ready_conf_cad <- df_new_with_cor_pre_taxo_cad |>
    make_other() |>
    no_other() |>
    prepare_plot_2()
  df_histogram_outside_ready_conf_cad <- df_peaks_samples_min_cad |>
    make_other() |>
    no_other() |>
    prepare_plot_2()

  df_histogram_ready_cad |>
    plot_histograms_taxo()
  # df_histogram_outside_ready_cad |>
  #   plot_histograms_taxo()

  df_histogram_ready_conf_cad |>
    plot_histograms_confident()
  #' keeping only distinct structures
  df_histogram_ready_conf_cad |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()
  df_histogram_outside_ready_conf_cad |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()

  df_histogram_ready_cad |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()
  df_histogram_outside_ready_cad |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()

  final_table_taxed_with_new_cor_cad <-
    df_new_with_cor_pre_taxo_cad |>
    prepare_hierarchy(detector = "cad")
  #' TODO check to use full table
  final_table_taxed_with_new_cor_ms_pos_cad <-
    df_new_with_cor_pre_taxo_cad |>
    prepare_hierarchy(detector = "ms")

  final_table_taxed_with_new_cor_conf_cad <-
    df_new_with_cor_pre_taxo_cad |>
    no_other() |>
    prepare_hierarchy(detector = "cad")
  final_table_taxed_with_new_cor_conf_ms_pos_cad <-
    df_new_with_cor_pre_taxo_cad |>
    no_other() |>
    prepare_hierarchy(detector = "ms")

  # samples_with_new_cor <-
  #   prepare_plot(dataframe = final_table_taxed_with_new_cor)
  #
  # absolute_with_new_cor <-
  #   plot_histograms(dataframe = samples_with_new_cor,
  #                   label = "CAD intensity of correlated peaks within CAD peak",
  #                   xlab = FALSE)
  #
  # absolute_with_new_cor

  plotly::plot_ly(
    data = final_table_taxed_with_new_cor_conf_cad,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_colors)

  plotly::plot_ly(
    data = final_table_taxed_with_new_cor_conf_ms_pos_cad,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_colors)

  index <- final_table_taxed_with_new_cor_cad |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))

  sunburst_grey_colors <- sunburst_colors
  sunburst_grey_colors[[which(index$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]

  plotly::plot_ly(
    data = final_table_taxed_with_new_cor_cad,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_grey_colors)

  plotly::plot_ly(
    data = final_table_taxed_with_new_cor_ms_pos_cad,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_grey_colors)
}

if (params$signal$detector$pda == TRUE) {
  log_debug(x = "loading compared peaks")
  compared_peaks_pda <-
    readr::read_delim(file = file.path(EXPORT_DIR, EXPORT_FILE_PDA))

  outside_peaks_pda <-
    readr::read_delim(file = file.path(EXPORT_DIR, EXPORT_FILE_PDA_2))

  log_debug(x = "joining compared peaks and candidates")
  df_peaks_samples_full_pda <- compared_peaks_pda |>
    dplyr::left_join(candidates_confident)

  df_peaks_samples_min_pda <- compared_peaks_pda |>
    dplyr::left_join(candidates_confident)

  log_debug(x = "temporary fix") #' TODO
  df_peaks_samples_full_pda <- df_peaks_samples_full_pda |>
    dplyr::mutate(
      id = sample,
      integral = peak_area,
      intensity = feature_area
    )
  df_peaks_samples_min_pda <- df_peaks_samples_min_pda |>
    dplyr::mutate(
      id = sample,
      integral = peak_area,
      intensity = feature_area,
      peak_rt_apex = rt,
      peak_area = 1
    )

  log_debug(x = "keeping peaks similarities above desired (pre-)threshold only")
  df_new_with_cor_pre_pda <- df_peaks_samples_full_pda |>
    dplyr::filter(comparison_score >= PEAK_SIMILARITY_PREFILTER) #' TODO check negative values

  log_debug(x = "keeping multiple features only if none was reported in the sepcies")
  df_new_with_cor_pre_taxo_pda <- df_new_with_cor_pre_pda |>
    dplyr::rowwise() |>
    dplyr::mutate(taxo = ifelse(
      test = grepl(
        pattern = best_candidate_organism,
        x = species
      ),
      yes = 1,
      no = 0
    )) |>
    dplyr::mutate(taxo = ifelse(
      test = is.na(taxo),
      yes = 0,
      no = taxo
    )) |>
    dplyr::mutate(taxo_2 = ifelse(
      test = best_candidate_organism %in% species,
      yes = 1,
      no = 0
    )) |>
    dplyr::group_by(sample, peak_id) |> #' TODO switch to ID if needed
    dplyr::mutate(sum = sum(taxo)) |>
    dplyr::mutate(sum_2 = sum(taxo_2)) |>
    dplyr::filter(taxo_2 == 1 |
      sum_2 == 0 |
      sum == 0) |> #' TODO decide if genus
    dplyr::ungroup()

  log_debug(x = "keeping peaks similarities with score above", PEAK_SIMILARITY)
  df_new_with_cor_pda <- df_new_with_cor_pre_taxo_pda |>
    dplyr::filter(comparison_score >= PEAK_SIMILARITY)

  log_debug(x = "plotting histograms")
  #' TODO harmonize 'others' among minor and major
  df_histogram_ready_pda <- df_new_with_cor_pre_taxo_pda |>
    make_other() |>
    prepare_plot_2()
  df_histogram_outside_ready_pda <- df_peaks_samples_min_pda |>
    make_other() |>
    prepare_plot_2()

  df_histogram_ready_conf_pda <- df_new_with_cor_pre_taxo_pda |>
    make_other() |>
    no_other() |>
    prepare_plot_2()
  df_histogram_outside_ready_conf_pda <- df_peaks_samples_min_pda |>
    make_other() |>
    no_other() |>
    prepare_plot_2()

  df_histogram_ready_pda |>
    plot_histograms_taxo()
  # df_histogram_outside_ready_pda |>
  #   plot_histograms_taxo()

  df_histogram_ready_conf_pda |>
    plot_histograms_confident()
  #' keeping only distinct structures
  df_histogram_ready_conf_pda |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()
  df_histogram_outside_ready_conf_pda |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()

  df_histogram_ready_pda |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()
  df_histogram_outside_ready_pda |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()

  final_table_taxed_with_new_cor_pda <-
    df_new_with_cor_pre_taxo_pda |>
    prepare_hierarchy(detector = "cad")
  #' TODO check to use full table
  final_table_taxed_with_new_cor_ms_pos_pda <-
    df_new_with_cor_pre_taxo_pda |>
    prepare_hierarchy(detector = "ms")

  final_table_taxed_with_new_cor_conf_pda <-
    df_new_with_cor_pre_taxo_pda |>
    no_other() |>
    prepare_hierarchy(detector = "cad")
  final_table_taxed_with_new_cor_conf_ms_pos_pda <-
    df_new_with_cor_pre_taxo_pda |>
    no_other() |>
    prepare_hierarchy(detector = "ms")

  # samples_with_new_cor <-
  #   prepare_plot(dataframe = final_table_taxed_with_new_cor)
  #
  # absolute_with_new_cor <-
  #   plot_histograms(dataframe = samples_with_new_cor,
  #                   label = "PDA intensity of correlated peaks within PDA peak",
  #                   xlab = FALSE)
  #
  # absolute_with_new_cor

  plotly::plot_ly(
    data = final_table_taxed_with_new_cor_conf_pda,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_colors)

  plotly::plot_ly(
    data = final_table_taxed_with_new_cor_conf_ms_pos_pda,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_colors)

  index <- final_table_taxed_with_new_cor_pda |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))

  sunburst_grey_colors <- sunburst_colors
  sunburst_grey_colors[[which(index$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]

  plotly::plot_ly(
    data = final_table_taxed_with_new_cor_pda,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_grey_colors)

  plotly::plot_ly(
    data = final_table_taxed_with_new_cor_ms_pos_pda,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_grey_colors)
}





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
