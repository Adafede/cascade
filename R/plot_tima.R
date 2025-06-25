#' Plot TIMA
#'
#' @export
#'
#' @include plot_histograms.R
#' @include prepare_hierarchy.R
#' @include prepare_plot.R
#' @include treemaps_progress.R
#'
#' @param tables Tables
#'
#' @return Pretty plots
#'
#' @examples NULL
plot_tima <- function(tables) {
  hierarchies <- tables |>
    purrr::map(
      prepare_hierarchy,
      type = "literature"
    )
  names(hierarchies) <- "comparison"

  prepared_plots <- hierarchies |>
    purrr::map(prepare_plot)

  histograms <- prepared_plots |>
    purrr::map(plot_histograms_litt, label = "")

  treemaps <-
    treemaps_progress_no_title(
      xs = names(hierarchies),
      hierarchies = hierarchies
    )

  sunbursts <-
    treemaps_progress_no_title(
      xs = names(hierarchies),
      type = "sunburst",
      hierarchies = hierarchies
    )

  return(list(
    histogram = histograms[[1]],
    treemap = treemaps[[1]],
    sunburst = sunbursts[[1]]
  ))
}
