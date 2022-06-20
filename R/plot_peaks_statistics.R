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
  
  plotly::plot_ly() |>
    plotly::add_pie(
      data = test_features,
      name = "Features",
      labels = ~ featuresPerPeak,
      values = ~ n,
      sort = FALSE,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 0)
    ) |>
    plotly::add_pie(
      data = test_structures,
      name = "Structures",
      labels = ~ structuresPerPeak,
      values = ~ n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 1)
    ) |>
    plotly::add_pie(
      data = test_classes,
      name = "Classes",
      labels = ~ chemicalClassesPerPeak,
      values = ~ n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 0, column = 2)
    ) |>
    plotly::add_pie(
      data = test_superclasses,
      name = "Superclasses",
      labels = ~ chemicalSuperclassesPerPeak,
      values = ~ n,
      type = "pie",
      textposition = "inside",
      domain = list(row = 1, column = 0)
    ) |>
    plotly::add_pie(
      data = test_pathways,
      name = "Pathways",
      labels = ~ chemicalPathwaysPerPeak,
      values = ~ n,
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
