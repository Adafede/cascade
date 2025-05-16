#' Plot TIMA histograms
#'
#' @export
#'
#' @include plot_histograms.R
#' @include prepare_hierarchy.R
#' @include prepare_plot.R
#'
#' @param tables Tables
#'
#' @return Pretty plots
#'
#' @examples NULL
plot_tima_histograms <- function(tables) {
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

  return(histograms$comparison)
}
