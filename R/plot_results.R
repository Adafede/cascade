#' Title
#'
#' @param detector
#'
#' @return
#' @export
#'
#' @examples
plot_results_1 <- function(detector = "cad") {
  list <- switch(detector,
    "bpi" = compared_peaks_list_bpi,
    "cad" = compared_peaks_list_cad,
    "pda" = compared_peaks_list_pda
  )

  log_debug(x = "preparing histograms")
  #' TODO harmonize 'others' among minor and major
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

  log_debug(x = "plotting histograms...")
  log_debug(x = "... taxo")
  histograms_taxo_maj <- df_histogram_maj |>
    plot_histograms_taxo()
  histograms_taxo_min <- df_histogram_min |>
    plot_histograms_taxo(level = "min")

  log_debug(x = "... confident features")
  histograms_conf_maj <- df_histogram_maj_conf |>
    plot_histograms_confident()
  histograms_conf_min <- df_histogram_min_conf |>
    plot_histograms_confident(level = "min")

  log_debug(x = "... unique structures")
  histograms_unique_maj <- df_histogram_maj |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()
  histograms_unique_min <- df_histogram_min |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident(level = "min")

  log_debug(x = "... confident structures")
  histograms_unique_conf_maj <- df_histogram_maj_conf |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident()
  histograms_unique_conf_min <- df_histogram_min_conf |>
    dplyr::distinct(peak_id, inchikey_2D, .keep_all = TRUE) |>
    plot_histograms_confident(level = "min")

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

#' Title
#'
#' @param detector
#'
#' @return
#' @export
#'
#' @examples
plot_results_2 <- function(detector = "cad") {
  list <- switch(detector,
    "bpi" = compared_peaks_list_bpi,
    "cad" = compared_peaks_list_cad,
    "pda" = compared_peaks_list_pda
  )

  temp_fix_duplicates <- function(df, colname = "peak_rt_apex") {
    df_new <- df |>
      dplyr::arrange(dplyr::desc(comparison_score)) |>
      dplyr::distinct(
        !!as.name(colname),
        inchikey_2D,
        best_candidate_1,
        best_candidate_2,
        best_candidate_3,
        .keep_all = TRUE
      )
    return(df_new)
  }
  temp_fix_posneg <- function(df) {
    df_new <- df |>
      dplyr::mutate(newrt = round(peak_rt_apex, 1)) |>
      dplyr::mutate(
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
      dplyr::mutate(
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
    return(df_new)
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

  log_debug(x = "preparing hierarchies...")
  log_debug(x = "... on everything")
  table_taxo_maj_cor_signal <- list$peaks_maj_precor_taxo_cor |>
    temp_fix_duplicates() |>
    prepare_hierarchy(detector = "cad")
  table_taxo_maj_cor_signal_2 <- list$peaks_maj_precor_taxo_cor |>
    temp_fix_posneg() |>
    temp_fix_duplicates(colname = "newrt") |>
    prepare_hierarchy(detector = "cad")


  #' TODO check to use full table
  table_taxo_maj_cor_ms <- list$peaks_maj_precor_taxo_cor |>
    prepare_hierarchy(detector = "ms")
  table_taxo_maj_cor_ms_2 <- table_taxo_maj_cor_ms |>
    dplyr::mutate(sample = gsub(
      pattern = "_pos",
      replacement = "",
      x = sample,
      ignore.case = TRUE
    )) |>
    dplyr::mutate(sample = gsub(
      pattern = "_neg",
      replacement = "",
      x = sample,
      ignore.case = TRUE
    )) |>
    dplyr::group_by(parents, ids, labels, species, sample) |>
    dplyr::summarise(values = sum(values))

  log_debug(x = "... on confident")
  table_taxo_maj_cor_conf_signal <-
    list$peaks_maj_precor_taxo_cor |>
    no_other() |>
    temp_fix_duplicates() |>
    prepare_hierarchy(detector = "cad")
  table_taxo_maj_cor_conf_signal_2 <-
    list$peaks_maj_precor_taxo_cor |>
    no_other() |>
    temp_fix_posneg() |>
    temp_fix_duplicates(colname = "newrt") |>
    prepare_hierarchy(detector = "cad")
  table_taxo_maj_cor_conf_ms <- list$peaks_maj_precor_taxo_cor |>
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

  #' TODO remove redundancy
  sunburst_conf_signal_based_pos <-
    table_taxo_maj_cor_conf_signal |>
    dplyr::filter(grepl(
      pattern = "_pos",
      x = sample,
      ignore.case = TRUE
    )) |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_colors)

  sunburst_conf_signal_based_neg <-
    table_taxo_maj_cor_conf_signal |>
    dplyr::filter(grepl(
      pattern = "_neg",
      x = sample,
      ignore.case = TRUE
    )) |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_colors)

  sunburst_conf_signal_based_duo <-
    table_taxo_maj_cor_conf_signal_2 |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_colors)

  sunburst_conf_ms_based_pos <- table_taxo_maj_cor_conf_ms |>
    dplyr::filter(grepl(
      pattern = "_pos",
      x = sample,
      ignore.case = TRUE
    )) |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_colors)

  sunburst_conf_ms_based_neg <- table_taxo_maj_cor_conf_ms |>
    dplyr::filter(grepl(
      pattern = "_neg",
      x = sample,
      ignore.case = TRUE
    )) |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_colors)

  index_signal_pos <- table_taxo_maj_cor_signal |>
    dplyr::filter(grepl(
      pattern = "_pos",
      x = sample,
      ignore.case = TRUE
    )) |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))
  index_signal_neg <- table_taxo_maj_cor_signal |>
    dplyr::filter(grepl(
      pattern = "_neg",
      x = sample,
      ignore.case = TRUE
    )) |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))
  index_signal_duo <- table_taxo_maj_cor_signal_2 |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))

  index_ms_pos <- table_taxo_maj_cor_ms |>
    dplyr::filter(grepl(
      pattern = "_pos",
      x = sample,
      ignore.case = TRUE
    )) |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))
  index_ms_neg <- table_taxo_maj_cor_ms |>
    dplyr::filter(grepl(
      pattern = "_neg",
      x = sample,
      ignore.case = TRUE
    )) |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))
  index_ms_duo <- table_taxo_maj_cor_ms_2 |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))

  sunburst_grey_colors_signal_pos <- sunburst_colors
  sunburst_grey_colors_signal_pos[[which(index_signal_pos$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]
  sunburst_grey_colors_ms_pos <- sunburst_colors
  sunburst_grey_colors_ms_pos[[which(index_ms_pos$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]
  sunburst_grey_colors_signal_neg <- sunburst_colors
  sunburst_grey_colors_signal_neg[[which(index_signal_neg$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]
  sunburst_grey_colors_ms_neg <- sunburst_colors
  sunburst_grey_colors_ms_neg[[which(index_ms_neg$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]
  sunburst_grey_colors_signal_duo <- sunburst_colors
  sunburst_grey_colors_signal_duo[[which(index_signal_duo$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]
  sunburst_grey_colors_ms_duo <- sunburst_colors
  sunburst_grey_colors_ms_duo[[which(index_ms_duo$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]

  sunburst_signal_based_pos <- table_taxo_maj_cor_signal |>
    dplyr::filter(grepl(
      pattern = "_pos",
      x = sample,
      ignore.case = TRUE
    )) |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_grey_colors_signal_pos)

  sunburst_signal_based_neg <- table_taxo_maj_cor_signal |>
    dplyr::filter(grepl(
      pattern = "_neg",
      x = sample,
      ignore.case = TRUE
    )) |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_grey_colors_signal_neg)

  sunburst_signal_based_duo <- table_taxo_maj_cor_signal_2 |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_grey_colors_signal_duo)

  sunburst_ms_based_pos <- table_taxo_maj_cor_ms |>
    dplyr::filter(grepl(
      pattern = "_pos",
      x = sample,
      ignore.case = TRUE
    )) |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_grey_colors_ms_pos)

  sunburst_ms_based_neg <- table_taxo_maj_cor_ms |>
    dplyr::filter(grepl(
      pattern = "_neg",
      x = sample,
      ignore.case = TRUE
    )) |>
    quick_sb() |>
    plotly::layout(colorway = sunburst_grey_colors_ms_neg)

  returned_list <- list(
    sunburst_signal_based_pos,
    sunburst_ms_based_pos,
    sunburst_conf_signal_based_pos,
    sunburst_conf_ms_based_pos,
    sunburst_signal_based_neg,
    sunburst_ms_based_neg,
    sunburst_conf_signal_based_neg,
    sunburst_conf_ms_based_neg,
    sunburst_signal_based_duo,
    sunburst_conf_signal_based_duo
  )
  names(returned_list) <- c(
    "sunburst_signal_based_pos",
    "sunburst_ms_based_pos",
    "sunburst_conf_signal_based_pos",
    "sunburst_conf_ms_based_pos",
    "sunburst_signal_based_neg",
    "sunburst_ms_based_neg",
    "sunburst_conf_signal_based_neg",
    "sunburst_conf_ms_based_neg",
    "sunburst_signal_based_duo",
    "sunburst_conf_signal_based_duo"
  )

  return(returned_list)
}
