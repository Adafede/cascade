#' TODO MOVE THAT PART AFTER
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

#' TODO Species step
#' 

#' TODO make confident
|>
  make_confident(score = CONFIDENCE_SCORE_MIN)


#' TODO rest

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