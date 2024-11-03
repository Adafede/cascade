#' Prepare plot
#'
#' @include colors.R
#'
#' @param dataframe Dataframe
#' @param organism Organism
#'
#' @return A dataframe prepared for plots
#'
#' @examples NULL
prepare_plot <- function(dataframe, organism = "species") {
  presamples <- dataframe |>
    tidytable::ungroup() |>
    tidytable::filter(parents != "" &
      !grepl(pattern = "-", x = parents)) |>
    tidytable::filter(parents != "" &
      !grepl(pattern = "-", x = parents)) |>
    tidytable::filter(!is.na(get(organism))) |>
    tidytable::mutate(species = get(organism)) |>
    tidytable::group_by(labels) |>
    tidytable::mutate(valuez = sum(values)) |>
    tidytable::ungroup() |>
    tidytable::arrange(desc(valuez)) |>
    tidytable::mutate(group = as.integer(factor(parents, levels = unique(parents)))) |>
    tidytable::group_by(group) |>
    tidytable::mutate(subgroup = as.integer(factor(x = ids, levels = unique(ids)))) |>
    tidytable::mutate(subgroup = tidytable::dense_rank(x = as.numeric(subgroup))) |>
    tidytable::rowwise() |>
    tidytable::group_by(sample) |>
    tidytable::mutate(tot = sum(values)) |>
    tidytable::ungroup()

  tempval <- presamples |>
    tidytable::filter(parents == "Other") |>
    tidytable::distinct(group) |>
    tidytable::pull()

  if (length(tempval) != 0) {
    samples_0 <- presamples |>
      tidytable::filter(group == tempval) |>
      tidytable::rowwise() |>
      tidytable::mutate(group = 666) |>
      tidytable::ungroup()

    samples_1 <- presamples |>
      tidytable::filter(group > tempval) |>
      tidytable::rowwise() |>
      tidytable::mutate(group = group - 1) |>
      tidytable::ungroup()

    samples_2 <- presamples |>
      tidytable::filter(group < tempval)

    samples <- rbind(samples_0, samples_1, samples_2)
  } else {
    samples <- presamples
  }

  samples <- samples |>
    tidytable::rowwise() |>
    tidytable::mutate(color = ifelse(
      test = grepl(
        pattern = "Other",
        x = parents,
        fixed = TRUE
      ),
      yes = microshades_grey[[1]][[subgroup]],
      no = microshades[[group]][[subgroup]]
    )) |>
    tidytable::mutate(relative = values / tot) |>
    tidytable::ungroup() |>
    tidytable::arrange(subgroup) |>
    tidytable::arrange(group) |>
    tidytable::mutate(idz = paste(group, subgroup))

  quickfix <- samples |>
    tidytable::filter(is.na(color)) |>
    tidytable::group_by(group, subgroup) |>
    tidytable::mutate(subsub = tidytable::cur_group_id()) |>
    tidytable::mutate(color = paste0("#00", sprintf("%04d", subsub))) |>
    tidytable::select(-subsub) |>
    tidytable::ungroup()

  samples <-
    rbind(samples |> tidytable::filter(!is.na(color)), quickfix)

  samples$ids <-
    forcats::fct_reorder2(
      .f = samples$ids,
      .x = samples$group,
      .y = samples$idz,
      .desc = FALSE
    )

  samples$color <-
    forcats::fct_reorder2(
      .f = samples$color,
      .x = samples$group,
      .y = samples$idz,
      .desc = FALSE
    )

  samples$sample <-
    forcats::fct_reorder2(
      .f = samples$sample,
      .x = samples$valuez,
      .y = samples$sample,
      .desc = TRUE
    )

  return(samples)
}

#' Prepare plot 2
#'
#' @include colors.R
#'
#' @param dataframe Dataframe
#'
#' @return A dataframe prepared for plots
#'
#' @examples NULL
prepare_plot_2 <- function(dataframe) {
  dataframe_prep <- dataframe |>
    tidytable::ungroup() |>
    tidytable::filter(!is.na(inchikey_2D)) |>
    tidytable::distinct(
      peak_id,
      inchikey_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      .keep_all = TRUE
    ) |>
    tidytable::group_by(peak_id) |>
    tidytable::mutate(peak_area = peak_area / max(tidytable::row_number())) |>
    tidytable::mutate(best_candidate_1 = ifelse(
      test = !is.na(best_candidate_1),
      yes = best_candidate_1,
      no = "Other"
    )) |>
    tidytable::mutate(
      best_candidate_1 = ifelse(
        test = best_candidate_1 != "notConfident",
        yes = best_candidate_1,
        no = "Other"
      )
    ) |>
    tidytable::mutate(
      best_candidate_1 = ifelse(
        test = best_candidate_1 != "notClassified",
        yes = best_candidate_1,
        no = "Other"
      )
    ) |>
    tidytable::mutate(
      best_candidate_1 = ifelse(
        test = best_candidate_1 != "notClassified",
        yes = best_candidate_1,
        no = "Other"
      )
    ) |>
    tidytable::mutate(best_candidate_2 = ifelse(
      test = !is.na(best_candidate_2),
      yes = best_candidate_2,
      no = "notAnnotated"
    )) |>
    tidytable::mutate(
      best_candidate_2 = ifelse(
        test = best_candidate_2 != "notConfident notConfident",
        yes = best_candidate_2,
        no = "notConfident"
      )
    ) |>
    tidytable::mutate(
      best_candidate_2 = ifelse(
        test = best_candidate_2 != "notClassified",
        yes = best_candidate_2,
        no = "notClassified"
      )
    ) |>
    tidytable::mutate(
      best_candidate_3 = ifelse(
        test = !is.na(best_candidate_3),
        yes = best_candidate_3,
        no = "notAnnotated notAnnotated"
      )
    ) |>
    tidytable::mutate(
      best_candidate_3 = ifelse(
        test = best_candidate_3 != "notConfident notConfident notConfident",
        yes = best_candidate_3,
        no = "notConfident notConfident"
      )
    ) |>
    tidytable::mutate(
      best_candidate_3 = ifelse(
        test = best_candidate_3 != "notClassified",
        yes = best_candidate_3,
        no = "notClassified notClassified"
      )
    ) |>
    tidytable::group_by(best_candidate_1) |>
    tidytable::mutate(sum = sum(unique(peak_area))) |>
    tidytable::ungroup() |>
    tidytable::arrange(tidytable::desc(sum)) |>
    tidytable::mutate(group = as.integer(factor(
      best_candidate_1,
      levels = unique(best_candidate_1)
    ))) |>
    tidytable::group_by(group) |>
    tidytable::mutate(subgroup = as.integer(factor(
      x = best_candidate_2, levels = unique(best_candidate_2)
    ))) |>
    tidytable::mutate(subgroup = tidytable::dense_rank(x = as.numeric(subgroup))) |>
    tidytable::rowwise() |>
    tidytable::mutate(color = ifelse(
      test = grepl(
        pattern = "Other",
        x = best_candidate_1,
        fixed = TRUE
      ),
      yes = microshades_grey[[1]][[subgroup]],
      no = microshades[[group]][[subgroup]]
    )) |>
    tidytable::mutate(name = paste(best_candidate_1, best_candidate_2, sep = " - ")) |>
    tidytable::select(
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
    tidytable::mutate(name_2 = ifelse(
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
    tidytable::mutate(name_2 = ifelse(
      test = is.na(name_2),
      yes = "Other",
      no = name_2
    )) |>
    tidytable::mutate(color_2 = ifelse(
      test = name_2 == "Species",
      yes = microshades[[9]][[1]],
      no = ifelse(
        test = name_2 == "Genus",
        yes = microshades[[9]][[2]],
        no = ifelse(
          test = name_2 == "Family",
          yes = microshades[[9]][[3]],
          no = ifelse(
            test = name_2 == "Kingdom",
            yes = microshades[[9]][[4]],
            no = microshades_grey[[1]][[2]]
          )
        )
      )
    )) |>
    tidytable::ungroup()

  if (nrow(dataframe_prep) > 0) {
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
  }

  return(dataframe_prep)
}
