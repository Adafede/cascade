df_histogram_test <- df_new_with_cor_08 |>
  dplyr::filter(id == "M_36_01") |>
  dplyr::distinct(
    peak_id,
    rt_apex,
    feature_id,
    integral,
    best_candidate_1,
    best_candidate_2,
    best_candidate_3
  ) |> 
  dplyr::mutate(
    best_candidate_1 = ifelse(
      test = !is.na(best_candidate_1),
      yes = best_candidate_1,
      no = "Other"
    )
  ) |>
  dplyr::mutate(
    best_candidate_2 = ifelse(
      test = !is.na(best_candidate_2),
      yes = best_candidate_2,
      no = "Other notClassified"
    )
  ) |>
  dplyr::mutate(
    best_candidate_3 = ifelse(
      test = !is.na(best_candidate_3),
      yes = best_candidate_3,
      no = "Other notClassified notClassified"
    )
  ) |> 
  dplyr::mutate(group = as.integer(factor(best_candidate_1, levels = unique(best_candidate_1)))) |>
  dplyr::group_by(group) |>
  dplyr::mutate(subgroup = as.integer(factor(x = best_candidate_2, levels = unique(best_candidate_2)))) |>
  dplyr::mutate(subgroup = dplyr::dense_rank(x = as.numeric(subgroup))) |>
  dplyr::rowwise() |>
  dplyr::mutate(color = ifelse(
    test = grepl(
      pattern = "Other",
      x = best_candidate_1,
      fixed = TRUE
    ),
    yes = grey_colors[[1]][subgroup],
    no = nice_colors[[group]][subgroup]
  ))

plotly::plot_ly(
  data = df_histogram_test,
  x = ~ rt_apex,
  y = ~ integral,
  color = ~ color,
  type = "bar"
)

df_histogram_test$color <-
  forcats::fct_reorder2(
    .f = df_histogram_test$color,
    .x = df_histogram_test$integral,
    .y = df_histogram_test$peak_id,
    .desc = TRUE
  )
ggplot2::ggplot(
  df_histogram_test,
  ggplot2::aes(
    x = rt_apex,
    y = integral,
    fill = best_candidate_2
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(
    values = levels(df_histogram_test$color) |>
      as.character(),
    guide = ggplot2::guide_legend(reverse = TRUE, ncol = 1)
  ) +
  ggplot2::labs(fill = "Chemical class") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.title = ggplot2::element_text(face = "bold"),
    # legend.position = "none",
    # axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
    axis.text.y = ggplot2::element_text(face = "italic")
  ) +
  ggplot2::ylab("Area")
