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
  df_histogram_min <- list$peaks_min |>
    make_other() |>
    prepare_plot_2()

  df_histogram_maj_conf <- list$peaks_maj_precor_taxo_cor |>
    make_other() |>
    no_other() |>
    prepare_plot_2()
  df_histogram_min_conf <- list$peaks_min |>
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

  log_debug(x = "preparing hierarchies...")
  log_debug(x = "... on everything")
  table_taxo_maj_cor_signal <- list$peaks_maj_precor_taxo_cor |>
    prepare_hierarchy(detector = "cad")
  #' TODO check to use full table
  table_taxo_maj_cor_ms <- list$peaks_maj_precor_taxo_cor |>
    prepare_hierarchy(detector = "ms")

  log_debug(x = "... on confident")
  table_taxo_maj_cor_conf_signal <-
    list$peaks_maj_precor_taxo_cor |>
    no_other() |>
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
  sunburst_conf_signal_based <- plotly::plot_ly(
    data = table_taxo_maj_cor_conf_signal,
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

  sunburst_conf_ms_based <- plotly::plot_ly(
    data = table_taxo_maj_cor_conf_ms,
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

  index_signal <- table_taxo_maj_cor_signal |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))
  index_ms <- table_taxo_maj_cor_ms |>
    dplyr::filter(parents == "") |>
    dplyr::arrange(desc(values))
  sunburst_grey_colors_signal <- sunburst_colors
  sunburst_grey_colors_signal[[which(index_signal$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]
  sunburst_grey_colors_ms <- sunburst_colors
  sunburst_grey_colors_ms[[which(index_ms$ids == "Other", arr.ind = TRUE)]] <-
    grey_colors[[1]][[5]]

  sunburst_signal_based <- plotly::plot_ly(
    data = table_taxo_maj_cor_signal,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_grey_colors_signal)

  sunburst_ms_based <- plotly::plot_ly(
    data = table_taxo_maj_cor_ms,
    ids = ~ids,
    labels = ~labels,
    parents = ~parents,
    values = ~values,
    maxdepth = 3,
    type = "sunburst",
    branchvalues = "total",
    textinfo = "label+percent value+percent parent+percent root"
  ) |>
    plotly::layout(colorway = sunburst_grey_colors_ms)

  returned_list <- list(
    sunburst_signal_based,
    sunburst_ms_based,
    sunburst_conf_signal_based,
    sunburst_conf_ms_based
  )
  names(returned_list) <- c(
    "sunburst_signal_based",
    "sunburst_ms_based",
    "sunburst_conf_signal_based",
    "sunburst_conf_ms_based"
  )

  return(returned_list)
}
