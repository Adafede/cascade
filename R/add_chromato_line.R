#' Add chromato line
#'
#' @param plot Plot
#' @param chromato Chromato
#' @param shift Shift
#' @param time Time
#' @param name Name
#' @param color Color
#' @param polarity Polarity
#'
#' @return A plot with added chromato line
#'
#' @examples NULL
add_chromato_line <- function(plot,
                              chromato,
                              shift = 0,
                              time,
                              name,
                              color,
                              polarity = "pos") {
  plot |>
    plotly::add_lines(
      data = chromato |>
        normalize_chromatograms_list(shift = shift, time = time),
      x = ~time,
      y = ~ switch(polarity,
        "pos" = intensity,
        "neg" = -1 * intensity
      ),
      name = name,
      line = list(width = 1, color = color)
    )
}
