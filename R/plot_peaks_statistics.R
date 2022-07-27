#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
plot_peaks_statistics <- function(df) {
  test_features <- df |>
    dplyr::group_by(featuresPerPeak) |>
    dplyr::count()
  test_mf <- df |>
    dplyr::group_by(molecularFormulasPerPeak) |>
    dplyr::count()
  test_structures <- df |>
    dplyr::group_by(structuresPerPeak) |>
    dplyr::count()
  test_classes <- df |>
    dplyr::group_by(chemicalClassesPerPeak) |>
    dplyr::count()
  test_superclasses <- df |>
    dplyr::group_by(chemicalSuperclassesPerPeak) |>
    dplyr::count()
  test_pathways <- df |>
    dplyr::group_by(chemicalPathwaysPerPeak) |>
    dplyr::count()

  test_features_cor <- df |>
    dplyr::group_by(correlatedFeaturesPerPeak) |>
    dplyr::count()
  test_mf_cor <- df |>
    dplyr::group_by(correlatedMolecularFormulasPerPeak) |>
    dplyr::count()
  test_structures_cor <- df |>
    dplyr::group_by(correlatedStructuresPerPeak) |>
    dplyr::count()
  test_classes_cor <- df |>
    dplyr::group_by(correlatedChemicalClassesPerPeak) |>
    dplyr::count()
  test_superclasses_cor <- df |>
    dplyr::group_by(correlatedChemicalSuperclassesPerPeak) |>
    dplyr::count()
  test_pathways_cor <- df |>
    dplyr::group_by(correlatedChemicalPathwaysPerPeak) |>
    dplyr::count()

  test_features_tax <- df |>
    dplyr::group_by(finalFeaturesPerPeak) |>
    dplyr::count()
  test_mf_tax <- df |>
    dplyr::group_by(finalMolecularFormulasPerPeak) |>
    dplyr::count()
  test_structures_tax <- df |>
    dplyr::group_by(finalStructuresPerPeak) |>
    dplyr::count()
  test_classes_tax <- df |>
    dplyr::group_by(finalChemicalClassesPerPeak) |>
    dplyr::count()
  test_superclasses_tax <- df |>
    dplyr::group_by(finalChemicalSuperclassesPerPeak) |>
    dplyr::count()
  test_pathways_tax <- df |>
    dplyr::group_by(finalChemicalPathwaysPerPeak) |>
    dplyr::count()

  test_full <- rbind(
    test_features,
    test_mf,
    test_structures,
    test_classes,
    test_superclasses,
    test_pathways
  ) |>
    dplyr::relocate(n, .before = 1) |>
    tidyr::pivot_longer(2:7) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::arrange(value) |>
    dplyr::group_by(name) |>
    dplyr::mutate(n = round(n / sum(n), 2))

  test_full_cor <- rbind(
    test_features_cor,
    test_mf_cor,
    test_structures_cor,
    test_classes_cor,
    test_superclasses_cor,
    test_pathways_cor
  ) |>
    dplyr::relocate(n, .before = 1) |>
    tidyr::pivot_longer(2:7) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::arrange(value) |>
    dplyr::group_by(name) |>
    dplyr::mutate(n = round(n / sum(n), 2))

  test_full_tax <- rbind(
    test_features_tax,
    test_mf_tax,
    test_structures_tax,
    test_classes_tax,
    test_superclasses_tax,
    test_pathways_tax
  ) |>
    dplyr::relocate(n, .before = 1) |>
    tidyr::pivot_longer(2:7) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::arrange(value) |>
    dplyr::group_by(name) |>
    dplyr::mutate(n = round(n / sum(n), 2))

  test_full <- within(
    test_full,
    name <- factor(
      name,
      levels = c(
        "featuresPerPeak",
        "molecularFormulasPerPeak",
        "structuresPerPeak",
        "chemicalClassesPerPeak",
        "chemicalSuperclassesPerPeak",
        "chemicalPathwaysPerPeak"
      )
    )
  )

  test_full_cor <- within(
    test_full_cor,
    name <- factor(
      name,
      levels = c(
        "correlatedFeaturesPerPeak",
        "correlatedMolecularFormulasPerPeak",
        "correlatedStructuresPerPeak",
        "correlatedChemicalClassesPerPeak",
        "correlatedChemicalSuperclassesPerPeak",
        "correlatedChemicalPathwaysPerPeak"
      )
    )
  )

  test_full_tax <- within(
    test_full_tax,
    name <- factor(
      name,
      levels = c(
        "finalFeaturesPerPeak",
        "finalMolecularFormulasPerPeak",
        "finalStructuresPerPeak",
        "finalChemicalClassesPerPeak",
        "finalChemicalSuperclassesPerPeak",
        "finalChemicalPathwaysPerPeak"
      )
    )
  )

  myDirtyColors <-
    c(
      "0" = "grey",
      "01" = "#33a02c",
      "02-05" = "#b2df8a",
      "06-10" = "#ff7f00",
      "10+" = "#fb9a99"
    )

  bars_1 <- ggplot2::ggplot(
    data = test_full,
    ggplot2::aes(
      x = name,
      y = n,
      fill = value,
      cumulative = TRUE
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = paste(value, ":", 100 * n, "%")),
      color = "white",
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::scale_colour_manual(
      values = myDirtyColors,
      aesthetics = c("colour", "fill")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      # axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

  bars_2 <- ggplot2::ggplot(
    data = test_full_cor,
    ggplot2::aes(
      x = name,
      y = n,
      fill = value,
      cumulative = TRUE
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = paste(value, ":", 100 * n, "%")),
      color = "white",
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::scale_colour_manual(
      values = myDirtyColors,
      aesthetics = c("colour", "fill")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      # axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

  bars_3 <- ggplot2::ggplot(
    data = test_full_tax,
    ggplot2::aes(
      x = name,
      y = n,
      fill = value,
      cumulative = TRUE
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = paste(value, ":", 100 * n, "%")),
      color = "white",
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::scale_colour_manual(
      values = myDirtyColors,
      aesthetics = c("colour", "fill")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      # axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

  circles_1 <- plotly::plot_ly(colors = myDirtyColors) |>
    plotly::add_pie(
      data = test_features_tax,
      name = "Features",
      labels = ~finalFeaturesPerPeak,
      values = ~n,
      sort = FALSE,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 0)
    ) |>
    plotly::add_pie(
      data = test_mf_tax,
      name = "Molecular Formulas",
      labels = ~finalMolecularFormulasPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 1)
    ) |>
    plotly::add_pie(
      data = test_structures_tax,
      name = "Structures",
      labels = ~finalStructuresPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 2)
    ) |>
    plotly::add_pie(
      data = test_classes_tax,
      name = "Classes",
      labels = ~finalChemicalClassesPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 1, column = 0)
    ) |>
    plotly::add_pie(
      data = test_superclasses_tax,
      name = "Superclasses",
      labels = ~finalChemicalSuperclassesPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 1, column = 1)
    ) |>
    plotly::add_pie(
      data = test_pathways_tax,
      name = "Pathways",
      labels = ~finalChemicalPathwaysPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 1, column = 2)
    ) |>
    plotly::layout(
      title = "Peak analysis \n Features > Structures > Chemical classes > \n Superclasses > Pathways",
      colorway = myDirtyColors,
      grid = list(rows = 2, columns = 3)
    )

  returned_list <- list(
    bars_1,
    # circles_1,
    bars_2,
    bars_3
  )
  names(returned_list) <- c(
    "bars_1",
    # circles_1,
    "bars_2",
    "bars_3"
  )
  return(returned_list)
}
