#' Plot peak detection
#'
#' @param df1 DF 1 containing chromatogram
#' @param df2 DF 2 containing peaks
#' @param fun Fun
#'
#' @return A plot with (non-)detected peaks
#'
#' @examples NULL
plot_peak_detection <- function(df1, df2, fun) {
  df1 |>
    plotly::plot_ly() |>
    plotly::add_trace(
      data = df1,
      x = ~time,
      y = ~intensity,
      type = "scatter",
      mode = "line",
      name = "signal",
      line = list(color = "1f78b4", width = 1)
    ) |>
    plotly::add_trace(
      data = df2,
      x = ~rt_apex,
      y = ~peak_max,
      yaxis = "y2",
      type = "scatter",
      marker = list(color = "ff7f00", symbol = "star"),
      name = "detected maximum",
      line = list(color = "1f78b4", width = 0)
    ) |>
    plotly::add_trace(
      data = df2,
      x = ~rt_min,
      y = ~ do.call(what = fun, args = list(rt_min)),
      yaxis = "y2",
      type = "scatter",
      marker = list(color = "33a02c", symbol = "triangle-right"),
      name = "detected minimum (start)",
      line = list(color = "1f78b4", width = 0)
    ) |>
    plotly::add_trace(
      data = df2,
      x = ~rt_max,
      y = ~ do.call(what = fun, args = list(rt_min)),
      yaxis = "y2",
      type = "scatter",
      marker = list(color = "33a02c", symbol = "triangle-left"),
      name = "detected minimum (end)",
      line = list(color = "1f78b4", width = 0)
    ) |>
    plotly::layout(
      yaxis = list(
        title = "Normalized Intensity",
        zeroline = TRUE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis2 = list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        overlaying = "y",
        side = "right"
      ),
      xaxis = list(
        title = "Time [min]",
        showticklabels = FALSE,
        showgrid = FALSE,
        overlaying = "y",
        side = "right"
      ),
      showlegend = FALSE
    )
}
