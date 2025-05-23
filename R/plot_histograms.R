#' Plot histograms
#'
#' @param dataframe Dataframe
#' @param chromatogram Chromatogram
#' @param label Label
#' @param y Y
#' @param xlab Xlab
#'
#' @return A plot of histograms
#'
#' @examples NULL
plot_histograms <-
  function(dataframe, chromatogram, label, y = "values", xlab = TRUE) {
    absolute <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = chromatogram,
        mapping = ggplot2::aes(x = time, y = intensity / max(intensity)),
        col = "black",
        size = 0.1
      ) +
      ggplot2::geom_bar(
        data = dataframe,
        mapping = ggplot2::aes(x = sample, y = get(y), fill = ids),
        color = "grey",
        stat = "identity",
        width = 1
      ) +
      ggplot2::scale_fill_manual(
        values = levels(dataframe$color) |>
          as.character(),
        guide = ggplot2::guide_legend()
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

#' Plot histograms confident
#'
#' @param dataframe Dataframe
#' @param chromatogram Chromatogram
#' @param level Level
#' @param time_min Time min
#' @param time_max Time max
#'
#' @return A plot of confident histograms
#'
#' @examples NULL
plot_histograms_confident <-
  function(dataframe, chromatogram, level = "max", time_min, time_max) {
    dataframe <- dataframe |>
      tidytable::group_by(peak_area) |>
      tidytable::mutate(n = max(tidytable::row_number())) |>
      tidytable::group_by(feature_area) |>
      tidytable::mutate(m = max(tidytable::row_number())) |>
      tidytable::ungroup()

    plot <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = chromatogram,
        mapping = ggplot2::aes(x = time, y = intensity / max(intensity)),
        col = "black",
        size = 0.1
      ) +
      {
        if (level == "max") {
          ggplot2::geom_bar(
            data = dataframe,
            mapping = ggplot2::aes(
              x = switch(level, "max" = peak_rt, "min" = feature_rt),
              y = switch(
                level,
                "max" = peak_area / max(peak_area * n),
                "min" = feature_area / max(feature_area * m)
              ),
              fill = name
            ),
            stat = "identity",
            width = 1
          ) # color = "grey",
        } else {
          ggplot2::geom_bar(
            data = dataframe,
            mapping = ggplot2::aes(
              x = switch(level, "max" = peak_rt, "min" = feature_rt),
              y = switch(
                level,
                "max" = peak_area / max(peak_area * n),
                "min" = feature_area / max(feature_area * m)
              ),
              fill = name
            ),
            stat = "identity",
            width = 1
          )
        }
      } +
      ggplot2::scale_fill_manual(
        values = levels(dataframe$color) |>
          as.character(),
        guide = ggplot2::guide_legend(ncol = 1)
      ) +
      ggplot2::labs(fill = "Chemical Pathway - Superclass") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.title = ggplot2::element_text(face = "bold"),
        # legend.title.align = 0.5,
        # legend.position = "bottom",
        # axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        # panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"),
        axis.text.y = ggplot2::element_text(face = "italic")
      ) +
      ggplot2::ylab("Intensity") +
      ggplot2::xlab("Retention time [min]") +
      ggplot2::xlim(max(time_min), min(time_max))

    return(plot)
  }

#' Plot histograms taxo
#'
#' @param dataframe Dataframe
#' @param chromatogram Chromatogram
#' @param level Level
#' @param mode Mode
#' @param time_min Time min
#' @param time_max Time max
#'
#' @return A plot of taxo histograms
#'
#' @examples NULL
plot_histograms_taxo <-
  function(
    dataframe,
    chromatogram,
    level = "max",
    mode = "pos",
    time_min,
    time_max
  ) {
    dataframe <- dataframe |>
      tidytable::group_by(peak_area) |>
      tidytable::mutate(n = max(tidytable::row_number())) |>
      tidytable::group_by(feature_area) |>
      tidytable::mutate(m = max(tidytable::row_number())) |>
      tidytable::ungroup()

    if (mode == "neg") {
      dataframe$peak_area <- -1 * dataframe$peak_area
      dataframe$feature_area <- -1 * dataframe$feature_area
      chromatogram$intensity <-
        -1 * chromatogram$intensity
    }

    plot <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = chromatogram,
        mapping = ggplot2::aes(x = time, y = intensity / max(abs(intensity))),
        col = "black",
        size = 0.1
      ) +
      {
        if (level == "max") {
          ggplot2::geom_bar(
            data = dataframe,
            mapping = ggplot2::aes(
              x = switch(level, "max" = peak_rt, "min" = feature_rt),
              y = switch(
                level,
                "max" = peak_area / max(abs(peak_area * n)),
                "min" = feature_area / max(abs(feature_area * m))
              ),
              fill = name_2
            ),
            stat = "identity",
            width = 1
          ) # color = "grey",
        } else {
          ggplot2::geom_bar(
            data = dataframe,
            mapping = ggplot2::aes(
              x = switch(level, "max" = peak_rt, "min" = feature_rt),
              y = switch(
                level,
                "max" = peak_area / max(abs(peak_area * n)),
                "min" = feature_area / max(abs(feature_area * m))
              ),
              fill = name_2
            ),
            stat = "identity",
            width = 1
          )
        }
      } +
      ggplot2::scale_fill_manual(
        values = levels(dataframe$color_2) |>
          as.character()
      ) +
      ggplot2::labs(fill = "Already reported in") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.title = ggplot2::element_text(face = "bold"),
        # legend.title.align = 0.5,
        # legend.position = "bottom",
        # axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        # panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"),
        axis.text.y = ggplot2::element_text(face = "italic")
      ) +
      ggplot2::ylab("Intensity") +
      ggplot2::xlab("Retention time [min]") +
      ggplot2::xlim(max(time_min), min(time_max))

    return(plot)
  }


#' Plot histograms litt
#'
#' @param dataframe Dataframe
#' @param label Label
#' @param y Y
#' @param xlab Xlab
#'
#' @return A plot of literature histograms
#'
#' @examples NULL
plot_histograms_litt <-
  function(dataframe, label, y = "values", xlab = TRUE) {
    absolute <- ggplot2::ggplot(
      dataframe,
      ggplot2::aes(x = sample, y = get(y), fill = ids)
    ) +
      ggplot2::geom_bar(
        stat = "identity",
        position = ggplot2::position_stack(reverse = TRUE)
      ) +
      ggplot2::scale_fill_manual(
        values = levels(dataframe$color) |>
          as.character(),
        guide = ggplot2::guide_legend(ncol = 1)
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
