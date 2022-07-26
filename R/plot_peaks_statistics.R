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

  test_full <- rbind(
    test_features,
    test_structures,
    test_classes,
    test_superclasses,
    test_pathways
  ) |>
    dplyr::relocate(n, .before = 1) |>
    tidyr::pivot_longer(2:6) |>
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
        "structuresPerPeak",
        "chemicalClassesPerPeak",
        "chemicalSuperclassesPerPeak",
        "chemicalPathwaysPerPeak"
      )
    )
  )

  ggplot2::ggplot(
    data = test_full,
    ggplot2::aes(
      x = name,
      y = n,
      fill = value,
      cumulative = TRUE
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = paste(value, ":", 100 * n, "%")),
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
}
