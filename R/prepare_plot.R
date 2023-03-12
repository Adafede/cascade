require(package = dplyr, quietly = TRUE)
require(package = forcats, quietly = TRUE)

#' Title
#'
#' @param dataframe
#' @param organism
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
    dplyr::mutate(species = get(organism)) |>
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
    dplyr::mutate(relative = values / tot) |>
    dplyr::ungroup()

  quickfix <- samples |>
    dplyr::filter(is.na(color)) |>
    dplyr::group_by(group, subgroup) |>
    dplyr::mutate(subsub = cur_group_id()) |>
    dplyr::mutate(color = paste0("#00", sprintf("%04d", subsub))) |>
    dplyr::select(-subsub) |>
    dplyr::ungroup()

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

#' Title
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
prepare_plot_2 <- function(dataframe) {
  dataframe_prep <- dataframe |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(inchikey_2D)) |>
    dplyr::distinct(
      peak_id,
      inchikey_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      .keep_all = TRUE
    ) |>
    dplyr::group_by(peak_id) |>
    dplyr::mutate(peak_area = peak_area / max(row_number())) |>
    dplyr::mutate(best_candidate_1 = ifelse(
      test = !is.na(best_candidate_1),
      yes = best_candidate_1,
      no = "Other"
    )) |>
    dplyr::mutate(
      best_candidate_1 = ifelse(
        test = best_candidate_1 != "notConfident",
        yes = best_candidate_1,
        no = "Other"
      )
    ) |>
    dplyr::mutate(
      best_candidate_1 = ifelse(
        test = best_candidate_1 != "notClassified",
        yes = best_candidate_1,
        no = "Other"
      )
    ) |>
    dplyr::mutate(
      best_candidate_1 = ifelse(
        test = best_candidate_1 != "notClassified",
        yes = best_candidate_1,
        no = "Other"
      )
    ) |>
    dplyr::mutate(best_candidate_2 = ifelse(
      test = !is.na(best_candidate_2),
      yes = best_candidate_2,
      no = "notAnnotated"
    )) |>
    dplyr::mutate(
      best_candidate_2 = ifelse(
        test = best_candidate_2 != "notConfident notConfident",
        yes = best_candidate_2,
        no = "notConfident"
      )
    ) |>
    dplyr::mutate(
      best_candidate_2 = ifelse(
        test = best_candidate_2 != "notClassified",
        yes = best_candidate_2,
        no = "notClassified"
      )
    ) |>
    dplyr::mutate(
      best_candidate_3 = ifelse(
        test = !is.na(best_candidate_3),
        yes = best_candidate_3,
        no = "notAnnotated notAnnotated"
      )
    ) |>
    dplyr::mutate(
      best_candidate_3 = ifelse(
        test = best_candidate_3 != "notConfident notConfident notConfident",
        yes = best_candidate_3,
        no = "notConfident notConfident"
      )
    ) |>
    dplyr::mutate(
      best_candidate_3 = ifelse(
        test = best_candidate_3 != "notClassified",
        yes = best_candidate_3,
        no = "notClassified notClassified"
      )
    ) |>
    dplyr::group_by(best_candidate_1) |>
    dplyr::mutate(sum = sum(unique(peak_area))) |>
    dplyr::ungroup() |>
    dplyr::arrange(desc(sum)) |>
    dplyr::mutate(group = as.integer(factor(
      best_candidate_1,
      levels = unique(best_candidate_1)
    ))) |>
    dplyr::group_by(group) |>
    dplyr::mutate(subgroup = as.integer(factor(
      x = best_candidate_2, levels = unique(best_candidate_2)
    ))) |>
    dplyr::mutate(subgroup = dplyr::dense_rank(x = as.numeric(subgroup))) |>
    dplyr::rowwise() |>
    dplyr::mutate(color = ifelse(
      test = grepl(
        pattern = "Other",
        x = best_candidate_1,
        fixed = TRUE
      ),
      yes = rev(grey_colors[[1]])[subgroup],
      no = rev(nice_colors[[group]])[subgroup]
    )) |>
    dplyr::mutate(name = paste(best_candidate_1, best_candidate_2, sep = " - ")) |>
    dplyr::select(
      sample,
      sample_organism = species,
      peak_id,
      # peak_rt_min,
      peak_rt = peak_rt_apex,
      # peak_rt_max,
      peak_area,
      feature_id,
      feature_area,
      feature_rt = rt,
      feature_mz = mz,
      comparison_score,
      smiles_2D,
      inchikey_2D,
      score_biological,
      score_chemical,
      score_final,
      nearest_organism = best_candidate_organism,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      # consensus_1,
      # consensus_2,
      # consensus_3,
      # consistency_1,
      # consistency_2,
      # consistency_3,
      # id,
      # integral,
      # intensity,
      # taxo,
      # sum,
      group,
      subgroup,
      color,
      name,
      sum
    ) |>
    dplyr::mutate(name_2 = ifelse(
      test = score_biological >= 0.9,
      yes = "Species",
      no = ifelse(
        test = score_biological >= 0.8,
        yes = "Genus",
        no = ifelse(
          test = score_biological >= 0.6,
          yes = "Family",
          no = ifelse(
            test = score_biological >= 0.2,
            yes = "Kingdom",
            no = "Other"
          )
        )
      )
    )) |>
    dplyr::mutate(name_2 = ifelse(
      test = is.na(name_2),
      yes = "Other",
      no = name_2
    )) |>
    dplyr::mutate(color_2 = ifelse(
      test = name_2 == "Species",
      yes = nice_colors[[9]][[5]],
      no = ifelse(
        test = name_2 == "Genus",
        yes = nice_colors[[9]][[4]],
        no = ifelse(
          test = name_2 == "Family",
          yes = nice_colors[[9]][[3]],
          no = ifelse(
            test = name_2 == "Kingdom",
            yes = nice_colors[[9]][[2]],
            no = grey_colors[[1]][[3]]
          )
        )
      )
    )) |>
    dplyr::ungroup()

  dataframe_prep$color <-
    forcats::fct_reorder2(
      .f = dataframe_prep$color,
      .x = dataframe_prep$sum,
      .y = dataframe_prep$group,
      .desc = TRUE
    )

  dataframe_prep$name <-
    forcats::fct_reorder2(
      .f = dataframe_prep$name,
      .x = dataframe_prep$sum,
      .y = dataframe_prep$group,
      .desc = TRUE
    )

  dataframe_prep$color_2 <-
    forcats::fct_reorder2(
      .f = dataframe_prep$color_2,
      .x = dataframe_prep$score_biological,
      .y = dataframe_prep$score_biological,
      .desc = TRUE
    )

  dataframe_prep$name_2 <-
    forcats::fct_reorder2(
      .f = dataframe_prep$name_2,
      .x = dataframe_prep$score_biological,
      .y = dataframe_prep$score_biological,
      .desc = TRUE
    )

  return(dataframe_prep)
}
