#' Plot chromatogram
#'
#' @noRd
#'
#' @param df Dataframe
#' @param text Text
#'
#' @return A plot of a chromatogram
#'
#' @examples NULL
# jarl-ignore unused_function: <reason>
plot_chromatogram <- function(df, text) {
  plotly::plot_ly(
    data = df,
    # |> tidytable::filter(grepl(pattern = "M", x = id)),
    x = ~time,
    y = ~intensity,
    color = ~id,
    colors = "Spectral",
    type = "scatter",
    mode = "lines",
    line = list(width = 0.5),
    legendgroup = ~id
  ) |>
    plotly::layout(
      annotations = list(
        x = 0.95,
        y = 0.95,
        xref = "paper",
        yref = "paper",
        text = text,
        showarrow = FALSE
      )
    )
}
