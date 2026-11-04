#' Plot TIMA annotation chemical class distributions
#'
#' Generates a stacked horizontal histogram, an interactive treemap, and an
#' interactive sunburst plot from one or more prepared TIMA annotation tables
#' returned by [prepare_tima_annotations()].
#'
#' @export
#'
#' @param tables A list of prepared annotation `data.table` objects (as
#'   returned by [prepare_tima_annotations()]).
#'
#' @return A named list with three elements:
#'   - `histogram`: A `ggplot` stacked bar chart.
#'   - `treemap`: A `plotly` interactive treemap widget.
#'   - `sunburst`: A `plotly` interactive sunburst widget.
#'
#' @concept EDAM:operation_3803 Natural product identification
#' @concept EDAM:operation_0337 Visualisation
#' @concept EDAM:operation_2942 Treemap visualisation
#' @concept EDAM:topic_3172 Metabolomics
#' @concept EDAM:topic_0092 Data visualisation
#'
#' @examples
#' tables <- prepare_tima_annotations(show_example = TRUE)
#' plot_tima(tables)
plot_tima <- function(tables) {
  hierarchies <- lapply(tables, prepare_hierarchy, type = "literature")
  names(hierarchies)[1L] <- "comparison"

  prepared_plots <- lapply(hierarchies, prepare_plot)
  histograms <- lapply(prepared_plots, plot_histograms_litt, label = "")

  treemaps <- treemaps_progress_no_title(
    xs = names(hierarchies),
    hierarchies = hierarchies
  )

  sunbursts <- treemaps_progress_no_title(
    xs = names(hierarchies),
    type = "sunburst",
    hierarchies = hierarchies
  )

  list(
    histogram = histograms[[1L]],
    treemap = treemaps[[1L]],
    sunburst = sunbursts[[1L]]
  )
}
