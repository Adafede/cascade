#' Plot results 1
#'
#' @include make_other.R
#' @include no_other.R
#' @include prepare_plot.R
#' @include plot_histograms.R
#'
#' @param list List
#' @param chromatogram Chromatogram
#' @param mode Mode
#' @param time_min Time min
#' @param time_max Time max
#'
#' @return A list of plots
#'
#' @examples NULL
plot_results_1 <- function(
  list,
  chromatogram,
  mode = "pos",
  time_min,
  time_max
) {
  message("preparing histograms")
  # TODO harmonize 'others' among minor and major
  df_histogram_maj <- list$peaks_maj_precor_taxo_cor |>
    make_other() |>
    prepare_plot_2()
  df_histogram_min <- list$peaks_min_precor_taxo_cor |>
    make_other() |>
    prepare_plot_2()

  df_histogram_maj_conf <- list$peaks_maj_precor_taxo_cor |>
    make_other() |>
    no_other() |>
    prepare_plot_2()
  df_histogram_min_conf <- list$peaks_min_precor_taxo_cor |>
    make_other() |>
    no_other() |>
    prepare_plot_2()

  message("plotting histograms...")
  message("... taxo")
  histograms_taxo_maj <- df_histogram_maj |>
    plot_histograms_taxo(
      chromatogram = chromatogram,
      mode = mode,
      time_min = time_min,
      time_max = time_max
    )
  histograms_taxo_min <- df_histogram_min |>
    plot_histograms_taxo(
      level = "min",
      chromatogram = chromatogram,
      mode = mode,
      time_min = time_min,
      time_max = time_max
    )

  message("... confident features")
  histograms_conf_maj <- df_histogram_maj_conf |>
    plot_histograms_confident(
      chromatogram = chromatogram,
      time_min = time_min,
      time_max = time_max
    )
  histograms_conf_min <- df_histogram_min_conf |>
    plot_histograms_confident(
      level = "min",
      chromatogram = chromatogram,
      time_min = time_min,
      time_max = time_max
    )

  message("... unique structures")
  histograms_unique_maj <- df_histogram_maj |>
    tidytable::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident(
      chromatogram = chromatogram,
      time_min = time_min,
      time_max = time_max
    )
  histograms_unique_min <- df_histogram_min |>
    tidytable::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident(
      level = "min",
      chromatogram = chromatogram,
      time_min = time_min,
      time_max = time_max
    )

  message("... confident structures")
  histograms_unique_conf_maj <- df_histogram_maj_conf |>
    tidytable::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident(
      chromatogram = chromatogram,
      time_min = time_min,
      time_max = time_max
    )
  histograms_unique_conf_min <- df_histogram_min_conf |>
    tidytable::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident(
      level = "min",
      chromatogram = chromatogram,
      time_min = time_min,
      time_max = time_max
    )

  returned_list <- list(
    histograms_taxo_maj,
    histograms_taxo_min,
    histograms_conf_maj,
    histograms_conf_min,
    histograms_unique_maj,
    histograms_unique_min,
    histograms_unique_conf_maj,
    histograms_unique_conf_min
  )
  names(returned_list) <- c(
    "histograms_taxo_maj",
    "histograms_taxo_min",
    "histograms_conf_maj",
    "histograms_conf_min",
    "histograms_unique_maj",
    "histograms_unique_min",
    "histograms_unique_conf_maj",
    "histograms_unique_conf_min"
  )

  return(returned_list)
}

#' Plot results 2
#'
#' @include prepare_hierarchy.R
#' @include no_other.R
#'
#' @param list List
#'
#' @return A list of plots
#'
#' @examples NULL
plot_results_2 <- function(list) {
  temp_fix_duplicates <- function(df, colname = "peak_rt_apex") {
    df |>
      tidytable::arrange(tidytable::desc(comparison_score)) |>
      tidytable::distinct(
        !!as.name(colname),
        inchikey_2D,
        best_candidate_1,
        best_candidate_2,
        best_candidate_3,
        .keep_all = TRUE
      )
  }

  temp_fix_posneg <- function(df) {
    df |>
      tidytable::mutate(newrt = round(peak_rt_apex, 1)) |>
      tidytable::mutate(
        sample = gsub(
          pattern = "_pos",
          replacement = "",
          x = sample,
          ignore.case = TRUE
        ),
        id = gsub(
          pattern = "_pos",
          replacement = "",
          x = id,
          ignore.case = TRUE
        )
      ) |>
      tidytable::mutate(
        sample = gsub(
          pattern = "_neg",
          replacement = "",
          x = sample,
          ignore.case = TRUE
        ),
        id = gsub(
          pattern = "_neg",
          replacement = "",
          x = id,
          ignore.case = TRUE
        )
      )
  }

  quick_sb <- function(df) {
    sb <- plotly::plot_ly(
      data = df,
      ids = ~ids,
      labels = ~labels,
      parents = ~parents,
      values = ~values,
      maxdepth = 3,
      type = "sunburst",
      branchvalues = "total",
      textinfo = "label+percent value+percent parent+percent root"
    )
    return(sb)
  }

  message("preparing hierarchies...")
  message("... on everything")
  # table_taxo_maj_cor_signal <- list$peaks_maj_precor_taxo_cor |>
  #   temp_fix_duplicates() |>
  #   prepare_hierarchy(detector = "cad")
  table_taxo_maj_cor_signal_2 <- list$peaks_maj_precor_taxo_cor |>
    temp_fix_posneg() |>
    temp_fix_duplicates(colname = "newrt") |>
    prepare_hierarchy(detector = "cad")

  ## TODO check to use full table
  table_taxo_maj_cor_ms <- list$peaks_maj_precor_taxo_cor |>
    prepare_hierarchy(detector = "ms")
  table_taxo_maj_cor_ms_2 <- table_taxo_maj_cor_ms |>
    tidytable::mutate(
      sample = gsub(
        pattern = "_pos",
        replacement = "",
        x = sample,
        ignore.case = TRUE
      )
    ) |>
    tidytable::mutate(
      sample = gsub(
        pattern = "_neg",
        replacement = "",
        x = sample,
        ignore.case = TRUE
      )
    ) |>
    tidytable::group_by(parents, ids, labels, species, sample) |>
    tidytable::summarize(values = sum(values))

  message("... on confident")
  # table_taxo_maj_cor_conf_signal <-
  #  list$peaks_maj_precor_taxo_cor |>
  #     no_other() |>
  #     temp_fix_duplicates() |>
  #     prepare_hierarchy(detector = "cad")
  table_taxo_maj_cor_conf_signal_1 <-
    list$peaks_maj_precor_taxo_cor |>
    no_other() |>
    temp_fix_posneg() |>
    temp_fix_duplicates(colname = "newrt")
  table_taxo_maj_cor_conf_signal_2 <-
    table_taxo_maj_cor_conf_signal_1 |>
    prepare_hierarchy(detector = "cad")
  # table_taxo_maj_cor_conf_ms <- list$peaks_maj_precor_taxo_cor |>
  #   no_other() |>
  #   prepare_hierarchy(detector = "ms")

  # samples_with_new_cor <-
  #   prepare_plot(dataframe = final_table_taxed_with_new_cor)
  #
  # absolute_with_new_cor <-
  #   plot_histograms(dataframe = samples_with_new_cor,
  #                   label = "CAD intensity of correlated peaks within CAD peak",
  #                   xlab = FALSE)
  #
  # absolute_with_new_cor

  ## TODO remove redundancy
  # sunburst_conf_signal_based_pos <-
  #   table_taxo_maj_cor_conf_signal |>
  #   tidytable::filter(grepl(
  #     pattern = "_pos",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   quick_sb() |>
  #   plotly::layout(colorway = microshades_colors)

  # sunburst_conf_signal_based_neg <-
  #   table_taxo_maj_cor_conf_signal |>
  #   tidytable::filter(grepl(
  #     pattern = "_neg",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   quick_sb() |>
  #   plotly::layout(colorway = microshades_colors)

  sunburst_conf_signal_based_duo <-
    table_taxo_maj_cor_conf_signal_2 |>
    quick_sb() |>
    plotly::layout(colorway = microshades_colors)

  # sunburst_conf_ms_based_pos <- table_taxo_maj_cor_conf_ms |>
  #   tidytable::filter(grepl(
  #     pattern = "_pos",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   quick_sb() |>
  #   plotly::layout(colorway = microshades_colors)

  # sunburst_conf_ms_based_neg <- table_taxo_maj_cor_conf_ms |>
  #   tidytable::filter(grepl(
  #     pattern = "_neg",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   quick_sb() |>
  #   plotly::layout(colorway = microshades_colors)

  # index_signal_pos <- table_taxo_maj_cor_signal |>
  #   tidytable::filter(grepl(
  #     pattern = "_pos",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   tidytable::filter(parents == "") |>
  #   tidytable::arrange(tidytable::desc(values))
  # index_signal_neg <- table_taxo_maj_cor_signal |>
  #   tidytable::filter(grepl(
  #     pattern = "_neg",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   tidytable::filter(parents == "") |>
  #   tidytable::arrange(tidytable::desc(values))

  index_signal_duo <- table_taxo_maj_cor_signal_2 |>
    tidytable::filter(parents == "") |>
    tidytable::arrange(tidytable::desc(values)) |>
    tidytable::distinct(ids, labels, .keep_all = TRUE)

  # index_ms_pos <- table_taxo_maj_cor_ms |>
  #   tidytable::filter(grepl(
  #     pattern = "_pos",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   tidytable::filter(parents == "") |>
  #   tidytable::arrange(tidytable::desc(values))
  # index_ms_neg <- table_taxo_maj_cor_ms |>
  #   tidytable::filter(grepl(
  #     pattern = "_neg",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   tidytable::filter(parents == "") |>
  #   tidytable::arrange(tidytable::desc(values))

  index_ms_duo <- table_taxo_maj_cor_ms_2 |>
    tidytable::filter(parents == "") |>
    tidytable::arrange(desc(values)) |>
    tidytable::distinct(ids, labels, .keep_all = TRUE)

  # sunburst_microshades_grey_signal_pos <- microshades_colors
  # index_s_pos <- which(index_signal_pos$ids == "Other", arr.ind = TRUE)
  # if (length(index_s_pos) > 0) {
  #   sunburst_microshades_grey_signal_pos[[index_s_pos]] <- microshades_grey[[1]][[5]]
  # }

  # sunburst_microshades_grey_ms_pos <- microshades_colors
  # index_m_pos <- which(index_ms_pos$ids == "Other", arr.ind = TRUE)
  # if (length(index_m_pos) > 0) {
  #   sunburst_microshades_grey_signal_pos[[index_m_pos]] <- microshades_grey[[1]][[5]]
  # }

  # sunburst_microshades_grey_signal_neg <- microshades_colors
  # index_s_neg <- which(index_signal_neg$ids == "Other", arr.ind = TRUE)
  # if (length(index_s_neg) > 0) {
  #   sunburst_microshades_grey_signal_pos[[index_s_neg]] <- microshades_grey[[1]][[5]]
  # }

  # sunburst_microshades_grey_ms_neg <- microshades_colors
  # index_m_neg <- which(index_ms_neg$ids == "Other", arr.ind = TRUE)
  # if (length(index_m_neg) > 0) {
  #   sunburst_microshades_grey_signal_pos[[index_m_neg]] <- microshades_grey[[1]][[5]]
  # }

  sunburst_microshades_grey_signal_duo <- microshades_colors
  index_s_duo <- which(index_signal_duo$ids == "Other", arr.ind = TRUE)
  if (length(index_s_duo) > 0) {
    sunburst_microshades_grey_signal_duo[[index_s_duo]] <- microshades_grey[[
      1
    ]][[5]]
  }

  # sunburst_microshades_grey_ms_duo <- microshades_colors
  # index_m_duo <- which(index_ms_duo$ids == "Other", arr.ind = TRUE)
  # if (length(index_m_duo) > 0) {
  #   sunburst_microshades_grey_signal_duo[[index_m_duo]] <- microshades_grey[[1]][[5]]
  # }

  # sunburst_signal_based_pos <- table_taxo_maj_cor_signal |>
  #   tidytable::filter(grepl(
  #     pattern = "_pos",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   quick_sb() |>
  #   plotly::layout(colorway = sunburst_microshades_grey_signal_pos)

  # sunburst_signal_based_neg <- table_taxo_maj_cor_signal |>
  #   tidytable::filter(grepl(
  #     pattern = "_neg",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   quick_sb() |>
  #   plotly::layout(colorway = sunburst_microshades_grey_signal_neg)

  sunburst_signal_based_duo <- table_taxo_maj_cor_signal_2 |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_microshades_grey_signal_duo)

  # sunburst_ms_based_pos <- table_taxo_maj_cor_ms |>
  #   tidytable::filter(grepl(
  #     pattern = "_pos",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   quick_sb() |>
  #   plotly::layout(colorway = sunburst_microshades_grey_ms_pos)

  # sunburst_ms_based_neg <- table_taxo_maj_cor_ms |>
  #   tidytable::filter(grepl(
  #     pattern = "_neg",
  #     x = sample,
  #     ignore.case = TRUE
  #   )) |>
  #   quick_sb() |>
  #   plotly::layout(colorway = sunburst_microshades_grey_ms_neg)

  returned_list <- list(
    # sunburst_signal_based_pos,
    # sunburst_ms_based_pos,
    # sunburst_conf_signal_based_pos,
    # sunburst_conf_ms_based_pos,
    # sunburst_signal_based_neg,
    # sunburst_ms_based_neg,
    # sunburst_conf_signal_based_neg,
    # sunburst_conf_ms_based_neg,
    sunburst_signal_based_duo,
    sunburst_conf_signal_based_duo
    # table_taxo_maj_cor_conf_signal_1
  )
  names(returned_list) <- c(
    # "sunburst_signal_based_pos",
    # "sunburst_ms_based_pos",
    # "sunburst_conf_signal_based_pos",
    # "sunburst_conf_ms_based_pos",
    # "sunburst_signal_based_neg",
    # "sunburst_ms_based_neg",
    # "sunburst_conf_signal_based_neg",
    # "sunburst_conf_ms_based_neg",
    "sunburst_signal_based_duo",
    "sunburst_conf_signal_based_duo"
    # "table_taxo_maj_cor_conf_signal_1"
  )

  return(returned_list)
}
