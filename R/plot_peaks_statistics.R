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
    scale_fill_viridis_c(option = "E") +
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
    scale_fill_viridis_c(option = "E") +
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

  circles_1 <- plotly::plot_ly() |>
    plotly::add_pie(
      data = test_features_cor,
      name = "Features",
      labels = ~correlatedFeaturesPerPeak,
      values = ~n,
      sort = FALSE,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 0)
    ) |>
    plotly::add_pie(
      data = test_mf_cor,
      name = "Molecular Formulas",
      labels = ~correlatedMolecularFormulasPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 1)
    ) |>
    plotly::add_pie(
      data = test_structures_cor,
      name = "Structures",
      labels = ~correlatedStructuresPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 2)
    ) |>
    plotly::add_pie(
      data = test_classes_cor,
      name = "Classes",
      labels = ~correlatedChemicalClassesPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 1, column = 0)
    ) |>
    plotly::add_pie(
      data = test_superclasses_cor,
      name = "Superclasses",
      labels = ~correlatedChemicalSuperclassesPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 1, column = 1)
    ) |>
    plotly::add_pie(
      data = test_pathways_cor,
      name = "Pathways",
      labels = ~correlatedChemicalPathwaysPerPeak,
      values = ~n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 1, column = 2)
    ) |>
    plotly::layout(
      title = "Peak analysis \n Features > Structures > Chemical classes > \n Superclasses > Pathways",
      colorway = viridis::cividis(max(
        nrow(test_features_cor),
        nrow(test_mf_cor),
        nrow(test_structures_cor),
        nrow(test_classes_cor),
        nrow(test_superclasses_cor),
        nrow(test_pathways_cor)
      )),
      grid = list(rows = 2, columns = 3)
    )

  returned_list <- list(
    bars_1,
    bars_2,
    circles_1
  )
  names(returned_list) <- c(
    "bars_1",
    "bars_2",
    "circles_1"
  )
  return(returned_list)
}
