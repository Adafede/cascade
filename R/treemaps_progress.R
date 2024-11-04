#' Treemaps progress
#'
#' @param xs XS
#' @param type Type
#' @param hierarchies Hierarchies
#'
#' @return A list of treemaps
#'
#' @examples NULL
treemaps_progress <- function(xs,
                              type = "treemap",
                              hierarchies) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = setNames(object = xs, nm = xs),
    FUN = function(x, hierarchies) {
      p()
      if (x != "special") {
        plotly::plot_ly(
          data = hierarchies[[x]],
          ids = ~ids,
          labels = ~labels,
          parents = ~parents,
          values = ~values,
          maxdepth = 3,
          type = type,
          branchvalues = "total",
          textinfo = "label+percent value+percent parent+percent root"
        ) |>
          plotly::layout(
            colorway = microshades_colors,
            title = paste(x, "(", nrow(
              tables[[x]] |> tidytable::distinct(structure)
            ), ")"),
            margin = list(t = 40)
          )
      } else {
        plotly::plot_ly() |>
          plotly::add_trace(
            data = hierarchies[[unique(hierarchies[[x]]$species)[1]]],
            ids = ~ids,
            labels = ~labels,
            parents = ~parents,
            values = ~values,
            maxdepth = 3,
            type = type,
            branchvalues = "total",
            textinfo = "label+percent value+percent parent+percent root",
            domain = list(row = 0, column = 0)
          ) |>
          plotly::add_trace(
            data = hierarchies[[unique(hierarchies[[x]]$species)[2]]],
            ids = ~ids,
            labels = ~labels,
            parents = ~parents,
            values = ~values,
            maxdepth = 3,
            type = type,
            branchvalues = "total",
            textinfo = "label+percent value+percent parent+percent root",
            domain = list(row = 0, column = 1)
          ) |>
          plotly::layout(
            title = paste(
              "Comparative analysis",
              "\n",
              unique(hierarchies[[x]]$species)[1],
              "(",
              nrow(tables[[unique(hierarchies[[x]]$species)[1]]] |>
                tidytable::distinct(structure)),
              ")",
              "                                 ",
              unique(hierarchies[[x]]$species)[2],
              "(",
              nrow(tables[[unique(hierarchies[[x]]$species)[2]]] |>
                tidytable::distinct(structure)),
              ")"
            ),
            grid = list(rows = 1, columns = 2),
            colorway = microshades_colors,
            margin = list(t = 40)
          )
      }
    },
    hierarchies = hierarchies
  )
}

#' Treemaps progress no title
#'
#' @param xs XS
#' @param type Type
#' @param hierarchies Hierarchies
#'
#' @return A list of treemaps with no title
#'
#' @examples NULL
treemaps_progress_no_title <- function(xs,
                                       type = "treemap",
                                       hierarchies) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = setNames(object = xs, nm = xs),
    FUN = function(x, hierarchies) {
      p()
      if (x != "special") {
        plotly::plot_ly(
          data = hierarchies[[x]],
          ids = ~ids,
          labels = ~labels,
          parents = ~parents,
          values = ~values,
          maxdepth = 3,
          type = type,
          branchvalues = "total",
          textinfo = "label+percent value+percent parent+percent root"
        ) |>
          plotly::layout(colorway = microshades_colors, margin = list(t = 40))
      } else {
        plotly::plot_ly() |>
          plotly::add_trace(
            data = hierarchies[[unique(hierarchies[[x]]$species)[1]]],
            ids = ~ids,
            labels = ~labels,
            parents = ~parents,
            values = ~values,
            maxdepth = 3,
            type = type,
            branchvalues = "total",
            textinfo = "label+percent value+percent parent+percent root",
            domain = list(row = 0, column = 0)
          ) |>
          plotly::add_trace(
            data = hierarchies[[unique(hierarchies[[x]]$species)[2]]],
            ids = ~ids,
            labels = ~labels,
            parents = ~parents,
            values = ~values,
            maxdepth = 3,
            type = type,
            branchvalues = "total",
            textinfo = "label+percent value+percent parent+percent root",
            domain = list(row = 0, column = 1)
          ) |>
          plotly::layout(
            grid = list(rows = 1, columns = 2),
            colorway = microshades_colors,
            margin = list(t = 40)
          )
      }
    },
    hierarchies = hierarchies
  )
}
