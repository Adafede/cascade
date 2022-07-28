#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
plot_peaks_statistics <- function(df) {
  accepted_variables <- c(
    "features",
    "molecular formulas",
    "structures",
    "chemical classes",
    "chemical superclasses",
    "chemical pathways"
  )
  myDirtyColors <-
    c(
      "0" = "grey",
      "01" = "#33a02c",
      "02-05" = "#b2df8a",
      "06-10" = "#ff7f00",
      "10+" = "#fb9a99"
    )

  inner_f <- function(variable) {
    leg <- paste("Number of", variable, "per peak")
    var <- variable |>
      gsub(pattern = " ", replacement = "") |>
      tolower()

    df_pretreated <- df |>
      tidyr::pivot_longer(cols = 14:19,
                          names_to = "names_3",
                          values_to = "3- peak shape + taxonomy + confidence filter") |>
      tidyr::pivot_longer(cols = 9:14,
                          names_to = "names_2",
                          values_to = "2- peak shape filter") |>
      tidyr::pivot_longer(cols = 3:8,
                          names_to = "names_1",
                          values_to = "1- no filter") |>
      tidyr::pivot_longer(cols = c(4, 6, 8), values_to = leg) |>
      dplyr::mutate_all(tolower) |>
      dplyr::rowwise() |>
      dplyr::filter(grepl(pattern = names_1,
                          x = names_2) &
                      grepl(pattern = names_1,
                            x = names_3))

    df_treated <- df_pretreated |>
      dplyr::filter(grepl(pattern = var,
                          x = names_1))

    alluvial <- ggplot2::ggplot(
      data = df_treated,
      mapping = ggplot2::aes(
        alluvium = peak_id,
        x = name,
        stratum = get(leg),
      )
    ) +
      geom_stratum(color = "black",
                   ggplot2::aes(fill = get(leg)),
                   decreasing = TRUE)  +
      geom_flow(
        ggplot2::aes(fill = get(leg),
                     colour = get(leg)),
        aes.bind = "alluvia",
        aes.flow = "backward",
        stat = ggplot2::after_stat("alluvium"),
        decreasing = TRUE
      ) +
      ggplot2::scale_fill_manual(values = myDirtyColors, name = leg) +
      ggplot2::scale_colour_manual(values = myDirtyColors, name = leg) +
      ggplot2::ylab("Number of peaks per category") +
      ggfittext::geom_fit_text(
        stat = "stratum",
        min.size = 0,
        fullheight = TRUE,
        reflow = TRUE,
        width = 1 / 4,
        ggplot2::aes(label = ggplot2::after_stat(stratum)),
        decreasing = TRUE
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank()
      )
  }
  alluvials_list <- lapply(X = accepted_variables, FUN = inner_f)

  return(alluvials_list)
}
