#' Title
#'
#' @param chromatograms
#' @param normalize_time
#' @param shift_cad
#' @param shift_pda
#' @param type
#'
#' @return
#' @export
#'
#' @examples
check_chromatograms <- function(chromatograms = c("bpi_pos", "cad_pos", "pda_pos"),
                                normalize_time = FALSE,
                                shift_cad = 0,
                                shift_pda = 0,
                                type = "improved") {
  stopifnot(
    "chromatograms must be in of 'bpi_pos', 'bpi_neg', 'cad_pos', 'cad_neg', 'pda_pos', 'pda_neg'" = chromatograms %in% c(
      "bpi_pos",
      "cad_pos",
      "pda_pos",
      "bpi_neg",
      "cad_neg",
      "pda_neg"
    )
  )
  stopifnot("type must be one of 'improved' or 'baselined'" = type %in% c("improved", "baselined"))
  plot <- plotly::plot_ly()
  if ("cad_pos" %in% chromatograms) {
    plot <- plot |>
      add_chromato_line(
        chromato = switch(type,
          "baselined" = chromatogram_cad_pos_baselined,
          "improved" =
            chromatogram_cad_pos_improved
        ),
        shift = shift_cad,
        time = normalize_time,
        name = "<b> CAD </b>",
        color = "#e31a1c"
      )
  }
  if ("cad_neg" %in% chromatograms) {
    plot <- plot |>
      add_chromato_line(
        chromato = switch(type,
          "baselined" = chromatogram_cad_neg_baselined,
          "improved" =
            chromatogram_cad_neg_improved
        ),
        shift = shift_cad,
        time = normalize_time,
        name = "<b> CAD (Neg)</b>",
        color = "#fb9a99",
        polarity = "neg"
      )
  }
  if ("pda_pos" %in% chromatograms) {
    plot <- plot |>
      add_chromato_line(
        chromato = switch(type,
          "baselined" = chromatogram_pda_pos_baselined,
          "improved" =
            chromatogram_pda_pos_improved
        ),
        shift = shift_pda,
        time = normalize_time,
        name = "<b> PDA </b>",
        color = "#b2df8a"
      )
  }
  if ("pda_neg" %in% chromatograms) {
    plot <- plot |>
      add_chromato_line(
        chromato = switch(type,
          "baselined" = chromatogram_pda_neg_baselined,
          "improved" =
            chromatogram_pda_neg_improved
        ),
        shift = shift_pda,
        time = normalize_time,
        name = "<b> PDA (Neg) </b>",
        color = "#33a02c",
        polarity = "neg"
      )
  }
  if ("bpi_pos" %in% chromatograms) {
    plot <- plot |>
      add_chromato_line(
        chromato = switch(type,
          "baselined" = chromatogram_bpi_pos_baselined,
          "improved" =
            chromatogram_bpi_pos_improved
        ),
        shift = 0,
        time = normalize_time,
        name = "<b> MS Pos </b>",
        color = "#a6cee3"
      )
  }
  if ("bpi_neg" %in% chromatograms) {
    plot <- plot |>
      add_chromato_line(
        chromato = switch(type,
          "baselined" = chromatogram_bpi_neg_baselined,
          "improved" =
            chromatogram_bpi_neg_improved
        ),
        shift = 0,
        time = normalize_time,
        name = "<b> MS Neg </b>",
        color = "#1f78b4",
        polarity = "neg"
      )
  }
  plot <- plot |>
    plotly::layout(
      xaxis = list(
        title = ifelse(
          normalize_time,
          "<b> Normalized Time </b>",
          "<b> Time [minutes] </b>"
        )
      ),
      yaxis = list(title = "<b> Normalized Intensity </b>")
    )
  plot
}
