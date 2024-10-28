#' Title
#'
#' @param plot
#' @param chromato
#' @param shift
#' @param time
#' @param name
#' @param color
#' @param polarity
#'
#' @return
#' @export
#'
#' @examples
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
