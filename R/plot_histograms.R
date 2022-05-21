require(package = ggplot2, quietly = TRUE)

#' Title
#'
#' @param dataframe
#' @param label
#' @param y
#' @param xlab
#'
#' @return
#' @export
#'
#' @examples
plot_histograms <-
  function(dataframe,
           label,
           y = "values",
           xlab = TRUE) {
    absolute <- ggplot2::ggplot(
      dataframe,
      ggplot2::aes(
        x = sample,
        y = get(y),
        fill = ids
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_bar(color = "grey", stat = "identity") +
      ggplot2::scale_fill_manual(
        values = levels(dataframe$color) |>
          as.character(),
        guide = ggplot2::guide_legend(reverse = TRUE, ncol = 1)
      ) +
      {
        if (xlab == TRUE) {
          ggplot2::scale_x_discrete(labels = levels(dataframe$sample))
        }
      } +
      ggplot2::labs(fill = "Chemical class") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.title = ggplot2::element_text(face = "bold"),
        # legend.position = "none",
        # axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        axis.text.y = ggplot2::element_text(face = "italic")
      ) +
      ggplot2::xlab(label) +
      ggplot2::ylab("Count") +
      ggplot2::coord_flip()

    return(absolute)
  }

#' Title
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
plot_histograms_confident <-
  function(dataframe) {
    plot <- ggplot2::ggplot(
      dataframe,
      ggplot2::aes(
        x = peak_rt,
        y = peak_area,
        fill = name
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_bar(color = "grey", stat = "identity") +
      ggplot2::scale_fill_manual(values = levels(dataframe$color) |>
        as.character()) +
      ggplot2::labs(fill = "Chemical Pathway - Superclass") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.title = ggplot2::element_text(face = "bold"),
        # legend.title.align = 0.5,
        # legend.position = "bottom",
        # axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = ggplot2::element_text(face = "italic")
      ) +
      ggplot2::ylab("Area") +
      ggplot2::xlab("Retention time [min]")

    return(plot)
  }

#' Title
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
plot_histograms_taxo <-
  function(dataframe) {
    plot <- ggplot2::ggplot(
      dataframe,
      ggplot2::aes(
        x = peak_rt,
        y = peak_area,
        fill = name_2
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_bar(color = "grey", stat = "identity") +
      ggplot2::scale_fill_manual(values = levels(dataframe$color_2) |>
        as.character()) +
      ggplot2::labs(fill = "Already reported in") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.title = ggplot2::element_text(face = "bold"),
        # legend.title.align = 0.5,
        # legend.position = "bottom",
        # axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.y = ggplot2::element_text(face = "italic")
      ) +
      ggplot2::ylab("Area") +
      ggplot2::xlab("Retention time [min]")

    return(plot)
  }
