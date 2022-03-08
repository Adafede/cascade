#' Title
#'
#' @param xs
#' @param type
#'
#' @return
#' @export
#'
#' @examples
treemaps_progress <- function(xs, type = "treemap") {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = setNames(
      object = xs,
      nm = xs
    ),
    FUN = function(x) {
      p()
      if (x != "special") {
        plotly::plot_ly(
          data = hierarchies[[x]] |>
            dplyr::filter(sample == x),
          ids = ~ids,
          labels = ~labels,
          parents = ~parents,
          values = ~values,
          maxdepth = 3,
          type = type,
          branchvalues = "total",
          textinfo = "label+value+percent parent+percent root"
        ) |>
          plotly::layout(
            colorway = sunburst_colors,
            title = paste(x, "(", nrow(
              tables[[x]] |> dplyr::distinct(structure)
            ), ")"),
            margin = list(t = 40)
          )
      } else {
        plotly::plot_ly() |>
          plotly::add_trace(
            data = hierarchies[[unique(hierarchies[[x]]$species)[1]]] |>
              dplyr::filter(sample == unique(hierarchies[[x]]$species)[1] &
                !is.na(labels)),
            ids = ~ids,
            labels = ~labels,
            parents = ~parents,
            values = ~values,
            maxdepth = 3,
            type = type,
            branchvalues = "total",
            textinfo = "label+value+percent parent+percent root",
            domain = list(row = 0, column = 0)
          ) |>
          plotly::add_trace(
            data = hierarchies[[unique(hierarchies[[x]]$species)[2]]] |>
              dplyr::filter(sample == unique(hierarchies[[x]]$species)[2] &
                !is.na(labels)),
            ids = ~ids,
            labels = ~labels,
            parents = ~parents,
            values = ~values,
            maxdepth = 3,
            type = type,
            branchvalues = "total",
            textinfo = "label+value+percent parent+percent root",
            domain = list(row = 0, column = 1)
          ) |>
          plotly::layout(
            title = paste(
              "Comparative analysis",
              "\n",
              unique(hierarchies[[x]]$species)[1],
              "(",
              nrow(tables[[unique(hierarchies[[x]]$species)[1]]] |> dplyr::distinct(structure)),
              ")",
              "                                 ",
              unique(hierarchies[[x]]$species)[2],
              "(",
              nrow(tables[[unique(hierarchies[[x]]$species)[2]]] |> dplyr::distinct(structure)),
              ")"
            ),
            grid = list(rows = 1, columns = 2),
            colorway = sunburst_colors,
            margin = list(t = 40)
          )
      }
    }
  )
}
