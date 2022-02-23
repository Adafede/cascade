require(package = dplyr, quietly = TRUE)
require(package = forcats, quietly = TRUE)

#' Title
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
prepare_plot <- function(dataframe, organism = "species") {
  presamples <- dataframe |>
    dplyr::ungroup() |>
    dplyr::filter(parents != "" &
      !grepl(pattern = "-", x = parents)) |>
    dplyr::filter(parents != "" &
      !grepl(pattern = "-", x = parents)) |>
    dplyr::filter(!is.na(get(organism))) |>
    dplyr::mutate(
      species = get(organism)
    ) |>
    dplyr::arrange(desc(values)) |>
    dplyr::mutate(group = as.integer(factor(parents, levels = unique(parents)))) |>
    dplyr::group_by(group) |>
    dplyr::mutate(subgroup = as.integer(factor(x = ids, levels = unique(ids)))) |>
    dplyr::mutate(subgroup = dplyr::dense_rank(x = as.numeric(subgroup))) |>
    dplyr::rowwise() |>
    dplyr::group_by(sample) |>
    dplyr::mutate(tot = sum(values)) |>
    dplyr::ungroup()

  tempval <- presamples |>
    dplyr::filter(parents == "Other") |>
    dplyr::distinct(group) |>
    dplyr::pull()

  if (length(tempval) != 0) {
    samples_0 <- presamples |>
      dplyr::filter(group == tempval) |>
      dplyr::rowwise() |>
      dplyr::mutate(group = 666) |>
      dplyr::ungroup()

    samples_1 <- presamples |>
      dplyr::filter(group > tempval) |>
      dplyr::rowwise() |>
      dplyr::mutate(group = group - 1) |>
      dplyr::ungroup()

    samples_2 <- presamples |>
      dplyr::filter(group < tempval)

    samples <- rbind(samples_0, samples_1, samples_2)
  } else {
    samples <- presamples
  }

  samples <- samples |>
    dplyr::rowwise() |>
    dplyr::mutate(color = ifelse(
      test = grepl(
        pattern = "Other",
        x = parents,
        fixed = TRUE
      ),
      yes = grey_colors[[1]][subgroup],
      no = nice_colors[[group]][subgroup]
    )) |>
    dplyr::mutate(relative = values / tot)

  quickfix <- samples |>
    dplyr::filter(is.na(color)) |>
    dplyr::group_by(group, subgroup) |>
    dplyr::mutate(subsub = cur_group_id()) |>
    dplyr::mutate(color = paste0("#00", sprintf("%04d", subsub))) |>
    dplyr::select(-subsub)

  samples <-
    rbind(samples |> dplyr::filter(!is.na(color)), quickfix)

  samples$ids <-
    forcats::fct_reorder2(
      .f = samples$ids,
      .x = samples$values,
      .y = samples$group,
      .desc = TRUE
    )

  samples$color <-
    forcats::fct_reorder2(
      .f = samples$color,
      .x = samples$values,
      .y = samples$group,
      .desc = TRUE
    )

  samples$sample <-
    forcats::fct_reorder2(
      .f = samples$sample,
      .x = samples$values,
      .y = samples$sample,
      .desc = FALSE
    )

  return(samples)
}
