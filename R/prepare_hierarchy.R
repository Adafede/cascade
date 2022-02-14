require(package = dplyr, quietly = TRUE)
require(package = splitstackshape, quietly = TRUE)

#' Title
#'
#' @param dataframe
#' @param type
#' @param detector
#' @param rescale
#'
#' @return
#' @export
#'
#' @examples
prepare_hierarchy <-
  function(dataframe,
           type = "analysis",
           detector = "ms",
           rescale = FALSE) {
    stopifnot("'type' must be either 'analysis' or 'literature'" = type %in% c("analysis", "literature"))

    if (type == "analysis") {
      df <- dataframe |>
        dplyr::mutate(best_candidate_1 = ifelse(
          test = !is.na(smiles_2D),
          yes = best_candidate_1,
          no = "Other"
        )) |>
        dplyr::mutate(best_candidate_2 = ifelse(
          test = !is.na(smiles_2D),
          yes = best_candidate_2,
          no = paste(best_candidate_1, "notAnnotated", sep = " ")
        )) |>
        dplyr::mutate(best_candidate_3 = ifelse(
          test = !is.na(smiles_2D),
          yes = best_candidate_3,
          no = paste(best_candidate_2, "notAnnotated", sep = " ")
        )) |>
        dplyr::mutate(
          best_candidate_2 = ifelse(
            test = best_candidate_2 != "notClassified",
            yes = best_candidate_2,
            no = paste(best_candidate_1, "notClassified", sep = " ")
          )
        ) |>
        dplyr::mutate(
          best_candidate_3 = ifelse(
            test = best_candidate_3 != "notClassified",
            yes = best_candidate_3,
            no = paste(best_candidate_2, "notClassified", sep = " ")
          )
        ) |>
        dplyr::mutate(
          best_candidate_1 = ifelse(
            test = best_candidate_3 != "notClassified notClassified notClassified",
            yes = best_candidate_1,
            no = "Other"
          )
        ) |>
        dplyr::mutate(
          best_candidate_2 = ifelse(
            test = best_candidate_3 != "notClassified notClassified notClassified",
            yes = best_candidate_2,
            no = "Other notClassified"
          )
        ) |>
        dplyr::mutate(
          best_candidate_3 = ifelse(
            test = best_candidate_3 != "notClassified notClassified notClassified",
            yes = best_candidate_3,
            no = "Other notClassified notClassified"
          )
        ) |>
        dplyr::mutate(chemical_pathway = best_candidate_1)

      df_notConfident <- df |>
        dplyr::filter(grepl(pattern = "notConfident", x = best_candidate_2)) |>
        dplyr::mutate(
          smiles_2D = NA,
          best_candidate_1 = "Other",
          best_candidate_2 = "Other notConfident",
          best_candidate_3 = "Other notConfident notConfident",
          inchikey_2D = NA
        ) |>
        dplyr::arrange(dplyr::desc(dplyr::across(dplyr::any_of(
          c("intensity", "comparison_score")
        )))) |>
        dplyr::distinct(id, peak_id, sample, .keep_all = TRUE)

      df_confident <- df |>
        dplyr::filter(!grepl(pattern = "notConfident", x = best_candidate_2))

      dataframe2 <- rbind(df_confident, df_notConfident) |>
        dplyr::group_by(chemical_pathway)
    } else {
      dataframe2 <- dataframe |>
        dplyr::mutate(
          best_candidate_1 = ifelse(
            test = !is.na(best_candidate_1),
            yes = best_candidate_1,
            no = "Other"
          ),
          best_candidate_2 = ifelse(
            test = !is.na(best_candidate_2),
            yes = best_candidate_2,
            no = paste(best_candidate_1, "notClassified", sep = " ")
          ),
          best_candidate_3 = ifelse(
            test = !is.na(best_candidate_3),
            yes = best_candidate_3,
            no = paste(best_candidate_2, "notClassified", sep = " ")
          )
        ) |>
        dplyr::mutate(chemical_pathway = ifelse(
          test = grepl(pattern = "not", x = best_candidate_1),
          yes = "Other",
          no = best_candidate_1
        ))
    }

    #' dirty temp fix because of bad npclassifier naming
    dataframe2 <- dataframe2 |>
      dplyr::rowwise() |>
      dplyr::mutate(safety = ifelse(
        test = best_candidate_3 == best_candidate_2,
        yes = "Y",
        no = "N"
      )) |>
      dplyr::mutate(best_candidate_3 = ifelse(
        test = safety == "Y",
        yes = paste(best_candidate_3, "sub"),
        no = best_candidate_3
      )) |>
      dplyr::group_by(chemical_pathway)

    parents <- dataframe2 |>
      dplyr::distinct(labels = best_candidate_1, ids = best_candidate_1) |>
      dplyr::mutate(parents = "") |>
      dplyr::distinct() |>
      dplyr::mutate(
        labels = gsub(
          pattern = "notAnnotated",
          replacement = "Other",
          x = labels,
          fixed = TRUE
        ),
        ids = gsub(
          pattern = "notAnnotated",
          replacement = "Other",
          x = ids,
          fixed = TRUE
        )
      ) |>
      dplyr::mutate(
        labels = gsub(
          pattern = "notClassified",
          replacement = "Other",
          x = labels,
          fixed = TRUE
        ),
        ids = gsub(
          pattern = "notClassified",
          replacement = "Other",
          x = ids,
          fixed = TRUE
        )
      ) |>
      dplyr::mutate(
        labels = gsub(
          pattern = "notConfident",
          replacement = "Other",
          x = labels,
          fixed = TRUE
        ),
        ids = gsub(
          pattern = "notConfident",
          replacement = "Other",
          x = ids,
          fixed = TRUE
        )
      ) |>
      dplyr::distinct()

    children_1 <- dataframe2 |>
      dplyr::distinct(best_candidate_1, best_candidate_2) |>
      dplyr::mutate(ids = paste(best_candidate_1, best_candidate_2, sep = "-")) |>
      dplyr::distinct(labels = best_candidate_2, ids) |>
      dplyr::mutate(parents = gsub(
        pattern = "-.*",
        replacement = "",
        x = ids
      )) |>
      dplyr::mutate_all(list(~ gsub(
        x = .x,
        pattern = "-NA$",
        replacement = "",
      ))) |>
      dplyr::filter(!is.na(labels)) |>
      dplyr::distinct()

    children_2 <- dataframe2 |>
      dplyr::distinct(best_candidate_2, best_candidate_3) |>
      dplyr::mutate(ids = paste(best_candidate_2, best_candidate_3, sep = "-")) |>
      dplyr::distinct(labels = best_candidate_3, ids, join = best_candidate_2) |>
      dplyr::mutate(
        join = gsub(
          pattern = "([a-zA-Z]|\\))(-)([a-zA-Z]|[0-9])(.*)",
          replacement = "\\1",
          x = join
        )
      ) |>
      dplyr::mutate_all(list(~ gsub(
        x = .x,
        pattern = "-NA$",
        replacement = "",
      ))) |>
      dplyr::mutate_all(list(~ gsub(
        x = .x,
        pattern = "^NA$",
        replacement = NA,
      ))) |>
      dplyr::full_join(children_1,
        by = c(
          "join" = "labels",
          "chemical_pathway" = "chemical_pathway"
        )
      ) |>
      dplyr::distinct(ids = ids.x, labels, parents = ids.y) |>
      dplyr::filter(!is.na(labels)) |>
      dplyr::distinct()

    genealogy <- rbind(parents, children_1, children_2) |>
      dplyr::ungroup() |>
      dplyr::distinct() |>
      dplyr::group_by(ids) |>
      dplyr::add_count() |> #' for ambiguous classes
      dplyr::filter(parents != "" |
        !is.na(ids) |
        !is.na(labels)) |>
      dplyr::ungroup()

    table <- dataframe2 |>
      dplyr::mutate(ids = paste0(best_candidate_2, "-", best_candidate_3)) |>
      dplyr::mutate(ids = gsub(
        x = ids,
        pattern = "-NA$",
        replacement = "",
      )) |>
      dplyr::mutate(ids = gsub(
        x = ids,
        pattern = "^NA$",
        replacement = NA,
      )) |>
      dplyr::full_join(genealogy) |>
      dplyr::mutate(best_candidate_3 = ifelse(
        test = is.na(best_candidate_3) & !is.na(best_candidate_2),
        no = best_candidate_3,
        yes = paste0(best_candidate_1, "-", best_candidate_2)
      ))

    if (type == "analysis") {
      table <- table |>
        dplyr::distinct(
          feature_id,
          smiles_2D,
          inchikey_2D,
          score_biological,
          score_chemical,
          score_final,
          best_candidate_organism,
          labels,
          ids,
          parents,
          sample,
          intensity,
          species,
          n
        )

      table_1 <- table |>
        dplyr::group_by(labels, sample) |>
        dplyr::add_count(name = "values") |>
        dplyr::select(
          chemical_pathway,
          parents,
          ids,
          labels,
          values,
          sample,
          intensity,
          species,
          n
        ) |>
        dplyr::distinct()

      table_1_1 <- table_1 |>
        dplyr::group_by(chemical_pathway, parents, sample) |>
        dplyr::mutate(values_2 = sum(values)) |>
        dplyr::ungroup()
    } else {
      table_1 <- table |>
        dplyr::group_by(labels, organism) |>
        dplyr::add_count(name = "values") |>
        dplyr::ungroup() |>
        dplyr::select(
          chemical_pathway,
          parents,
          ids,
          labels,
          values,
          organism,
          sample,
          n
        ) |>
        dplyr::distinct()

      table_1_1 <- table_1 |>
        dplyr::group_by(chemical_pathway, parents, organism) |>
        dplyr::mutate(values_2 = sum(values)) |>
        dplyr::ungroup()
    }

    top_parents_table <-
      dplyr::left_join(
        table_1_1 |> distinct(chemical_pathway, parents, ids, labels, values_2, n),
        table_1_1 |> distinct(chemical_pathway, parents, ids, labels, values_2, n),
        by = c("ids" = "parents")
      ) |>
      dplyr::filter(!grepl(pattern = "-", x = parents)) |>
      dplyr::distinct(
        chemical_pathway = chemical_pathway.x,
        parents,
        ids,
        labels = labels.x,
        values_2 = values_2.y,
        n.x
      ) |>
      dplyr::group_by(chemical_pathway, parents) |>
      dplyr::mutate(values_3 = sum(!is.na(values_2)) / as.numeric(n.x)) |>
      dplyr::filter(parents != "") |>
      dplyr::distinct(chemical_pathway, parents, values_3) |>
      dplyr::ungroup() |>
      dplyr::top_n(n = length(nice_colors), wt = values_3) |>
      dplyr::arrange(desc(values_3)) |>
      dplyr::select(-values_3) |>
      dplyr::mutate(ids = parents, labels = parents) |>
      dplyr::mutate(parents = "", new_labels = labels) |>
      dplyr::distinct()

    if (nrow(top_parents_table > length(nice_colors))) {
      top_parents_table <- top_parents_table |>
        head(length(nice_colors)) #' in case of equal numbers among classes
    }

    top_parents <- top_parents_table$labels

    table_2 <- table_1_1 |>
      dplyr::group_by(chemical_pathway, parents) |>
      dplyr::summarize(sum = sum(values))

    table_3 <-
      dplyr::left_join(table_2, children_1, by = c("parents" = "ids")) |>
      dplyr::group_by(chemical_pathway = chemical_pathway.y, parents.y) |>
      dplyr::mutate(n = sum(sum)) |>
      dplyr::mutate(n = if_else(
        condition = parents.y == "Other",
        true = n + 666,
        false = n + 0
      )) |>
      dplyr::select(parents = parents.y, ids = parents, labels, sum, n) |>
      dplyr::ungroup()

    table_3_1 <- table_3 |>
      dplyr::distinct(n, .keep_all = TRUE) |>
      dplyr::slice_max(n, n = 4, with_ties = FALSE) |>
      dplyr::select(
        chemical_pathway,
        parents,
        ids,
        labels
      )

    top_medium_table <- table_3 |>
      dplyr::filter(parents %in% top_parents) |>
      dplyr::group_by(parents) |>
      dplyr::distinct(sum, .keep_all = TRUE) |>
      dplyr::slice_max(sum, n = 4, with_ties = FALSE) |>
      dplyr::select(
        chemical_pathway,
        parents,
        ids,
        labels
      ) |>
      dplyr::mutate(new_labels = labels) |>
      dplyr::distinct()

    top_medium <- top_medium_table$ids

    low_medium_table <-
      dplyr::anti_join(table_3, dplyr::bind_rows(table_3_1, top_medium_table)) |>
      dplyr::select(
        chemical_pathway,
        parents,
        ids,
        labels
      ) |>
      dplyr::filter(!is.na(parents)) |>
      dplyr::distinct()

    low_medium_table_1_a <- low_medium_table |>
      dplyr::filter(parents %in% top_parents) |>
      dplyr::mutate(new_labels = paste(parents, "Other")) |>
      dplyr::mutate(ids = paste(parents, new_labels, sep = "-")) |>
      dplyr::mutate(labels = new_labels) |>
      dplyr::distinct()

    low_medium_table_1_b <- low_medium_table |>
      dplyr::filter(parents %in% top_parents) |>
      dplyr::mutate(new_labels = paste(parents, "Other")) |>
      dplyr::mutate(
        parents = paste(parents, new_labels, sep = "-"),
        ids = paste(new_labels, labels, sep = "-")
      ) |>
      dplyr::mutate(new_labels = labels) |>
      dplyr::distinct()

    low_medium_table_1 <-
      dplyr::bind_rows(low_medium_table_1_a, low_medium_table_1_b) |>
      dplyr::distinct()

    low_medium_table_2 <- low_medium_table |>
      dplyr::filter(!parents %in% top_parents) |>
      dplyr::mutate(new_labels = paste("Other", parents)) |>
      dplyr::mutate(parents = "Other") |>
      dplyr::mutate(ids = paste(parents, new_labels, sep = "-")) |>
      dplyr::distinct()

    genealogy_new_med <-
      dplyr::bind_rows(
        top_parents_table,
        top_medium_table,
        low_medium_table_2,
        low_medium_table_1
      ) |>
      dplyr::distinct()

    table_next <- genealogy |>
      dplyr::filter(!labels %in% genealogy_new_med$labels) |>
      dplyr::filter(parents %in% top_medium_table$ids) |>
      dplyr::filter(parents != "") |>
      dplyr::mutate(new_labels = paste(
        gsub(
          pattern = "^([^-]*)-(.*)",
          replacement = "",
          x = parents
        ),
        labels
      ) |>
        trimws()) |>
      dplyr::mutate(ids = paste(
        gsub(
          pattern = "^([^-]*)-(.*)",
          replacement = "\\2",
          x = parents
        ),
        new_labels,
        sep = "-"
      )) |>
      dplyr::mutate(new_labels = if_else(
        condition = grepl(pattern = "^Other", x = parents),
        true = labels,
        false = new_labels
      ))

    genealogy_new_med_2 <-
      dplyr::bind_rows(genealogy_new_med, table_next) |>
      dplyr::ungroup() |>
      dplyr::distinct()

    table_children <- genealogy |>
      dplyr::filter(!labels %in% genealogy_new_med$ids) |>
      dplyr::select(labels) |>
      dplyr::left_join(genealogy_new_med_2) |>
      dplyr::filter(grepl(pattern = "-", x = parents)) |>
      dplyr::left_join(genealogy, by = c("parents" = "ids")) |>
      dplyr::select(
        chemical_pathway = chemical_pathway.x,
        parents,
        ids,
        labels = labels.x,
        new_labels
      )

    table_other <- genealogy_new_med_2 |>
      dplyr::filter(grepl(pattern = " Other-", x = ids)) |>
      dplyr::distinct(chemical_pathway, ids = parents) |>
      dplyr::mutate(parents = (gsub(
        pattern = "-.*",
        replacement = "\\1",
        x = ids
      ))) |>
      dplyr::mutate(labels = (gsub(
        pattern = ".*-",
        replacement = "\\1",
        x = ids
      )), new_labels = labels)

    genealogy_new_med_3 <-
      dplyr::bind_rows(
        genealogy_new_med_2,
        table_children,
        table_other
      ) |>
      dplyr::distinct()

    missing_children_1 <- genealogy_new_med_3 |>
      dplyr::filter(grepl(pattern = " Other$", x = parents)) |>
      dplyr::distinct(chemical_pathway,
        parents = ids,
        best_candidate_2 = labels
      ) |>
      dplyr::left_join(dataframe2 |>
        dplyr::distinct(best_candidate_3, best_candidate_2)) |>
      dplyr::mutate(
        ids = paste(best_candidate_2,
          best_candidate_3,
          sep = "-"
        ),
        labels = best_candidate_3,
        new_labels = best_candidate_3
      ) |>
      dplyr::select(chemical_pathway, parents, ids, labels, new_labels) |>
      dplyr::mutate_all(list(~ gsub(
        x = .x,
        pattern = "-NA$",
        replacement = "",
      )))

    missing_children_2 <- genealogy_new_med_3 |>
      dplyr::filter(grepl(pattern = "^Other", x = parents)) |>
      dplyr::distinct(chemical_pathway,
        parents = ids,
        best_candidate_2 = labels
      ) |>
      dplyr::left_join(dataframe2 |>
        dplyr::distinct(best_candidate_3, best_candidate_2)) |>
      dplyr::mutate(
        ids = paste(
          gsub(
            pattern = ".*-",
            replacement = "\\1",
            x = parents
          ),
          best_candidate_3,
          sep = "-"
        ),
        labels = best_candidate_3,
        new_labels = best_candidate_3
      ) |>
      dplyr::select(chemical_pathway, parents, ids, labels, new_labels) |>
      dplyr::distinct() |>
      dplyr::mutate_all(list(~ gsub(
        x = .x,
        pattern = "-NA$",
        replacement = "",
      ))) |>
      dplyr::filter(!is.na(labels))

    missing_children <-
      rbind(missing_children_1, missing_children_2)

    genealogy_new_med_4 <-
      dplyr::bind_rows(genealogy_new_med_3, missing_children) |>
      dplyr::distinct(chemical_pathway, parents, ids, labels, new_labels)

    table_new <- dataframe2 |>
      dplyr::filter(!is.na(species)) |>
      dplyr::full_join(
        genealogy_new_med_4,
        by = c(
          "best_candidate_3" = "labels",
          "chemical_pathway" = "chemical_pathway"
        )
      )

    if (type == "analysis") {
      table_new <- table_new |>
        dplyr::rowwise() |>
        dplyr::filter(grepl(pattern = id, x = sample) |
          #' Comment these if
          is.na(id)) |> #' using wrong feature table
        dplyr::ungroup() |>
        dplyr::mutate(intensity = switch(detector,
          "ms" = intensity,
          "cad" = integral
        )) |>
        dplyr::distinct(
          feature_id,
          smiles_2D,
          inchikey_2D,
          score_biological,
          score_chemical,
          score_final,
          best_candidate_organism,
          labels = new_labels,
          ids,
          parents,
          sample,
          intensity,
          species
        )

      table_1_new <- table_new |>
        dplyr::group_by(parents, ids, sample, intensity) |>
        dplyr::add_count(name = "values") |>
        dplyr::select(parents, ids, labels, values, sample, intensity, species) |>
        dplyr::distinct() |>
        dplyr::filter(!is.na(ids))
    } else {
      table_1_new <- table_new |>
        dplyr::mutate(labels = best_candidate_3) |>
        dplyr::group_by(parents, ids, sample) |>
        dplyr::add_count(name = "values") |>
        dplyr::select(parents, ids, labels, values, sample, organism, species) |>
        dplyr::distinct() |>
        dplyr::filter(!is.na(ids))
    }

    additional_row <- table_1_new |>
      dplyr::filter(labels == "Other other") |>
      dplyr::mutate(
        parents = "Other-Other other",
        ids = "Other other - Other other other",
        labels = "Other other other"
      )

    table_1_1_new <- rbind(table_1_new, additional_row) |>
      dplyr::ungroup()

    final_table_3_1 <- table_1_1_new |>
      dplyr::filter(ids %in% low_medium_table_1_a$ids) |>
      dplyr::distinct(parents, ids, labels)

    final_table_4_1 <- table_1_1_new |>
      dplyr::filter(grepl(pattern = "-", x = parents) &
        parents %in% missing_children_1$parents) |>
      dplyr::group_by(dplyr::across(dplyr::any_of(
        c("parents", "ids", "labels", "sample", "species")
      ))) |>
      dplyr::summarise(values = sum(switch(type,
        "analysis" = intensity,
        "literature" = values
      ))) |>
      dplyr::ungroup() |>
      dplyr::filter(!is.na(labels)) |>
      dplyr::mutate(join = gsub(
        pattern = "-.*",
        replacement = "",
        x = parents
      ))

    final_table_4 <- final_table_4_1 |>
      dplyr::select(-join)

    final_table_3_2 <- final_table_3_1 |>
      dplyr::left_join(final_table_4_1, by = c("labels" = "join")) |>
      dplyr::select(
        parents = ids.x,
        ids = parents.y,
        labels = ids.y,
        sample,
        species,
        values
      ) |>
      dplyr::mutate(
        labels = gsub(
          pattern = "([a-zA-Z]|\\))(-)([a-zA-Z]|[0-9])(.*)",
          replacement = "\\1",
          x = labels
        )
      ) |>
      dplyr::group_by(dplyr::across(dplyr::any_of(
        c("parents", "ids", "labels", "sample", "species")
      ))) |>
      dplyr::summarise(values = sum(values)) |>
      dplyr::ungroup()

    final_table_3_3 <- table_1_1_new |>
      dplyr::filter(!is.na(species)) |>
      dplyr::filter(grepl(
        pattern = "-",
        x = parents
      ) &
        !parents %in% missing_children_1$parents) |>
      dplyr::group_by(dplyr::across(dplyr::any_of(
        c("parents", "ids", "labels", "sample", "species")
      ))) |>
      dplyr::summarise(values = sum(switch(type,
        "analysis" = intensity,
        "literature" = values
      ))) |>
      dplyr::ungroup()

    final_table_3 <-
      dplyr::bind_rows(
        final_table_3_3,
        final_table_3_2
      ) |>
      dplyr::distinct() |>
      dplyr::arrange(parents, ids)

    final_table_2 <- table_1_1_new |>
      dplyr::select(-any_of(c(
        "values", "sample", "intensity", "species", "organism"
      ))) |>
      dplyr::left_join(final_table_3 |>
        dplyr::select(any_of(
          c("values", "sample", "parents", "species")
        )),
      by = c("ids" = "parents")
      ) |>
      dplyr::group_by(dplyr::across(dplyr::any_of(
        c("parents", "ids", "labels", "sample", "species")
      ))) |>
      dplyr::summarise(values = sum(values)) |>
      dplyr::ungroup() |>
      dplyr::filter(!is.na(values))

    final_table_1 <- table_1_1_new |>
      dplyr::select(-any_of(c(
        "values", "sample", "intensity", "species", "organism"
      ))) |>
      dplyr::left_join(final_table_2 |>
        dplyr::select(any_of(
          c("values", "sample", "parents", "species")
        )),
      by = c("ids" = "parents")
      ) |>
      dplyr::group_by(dplyr::across(dplyr::any_of(
        c("parents", "ids", "labels", "sample", "species")
      ))) |>
      dplyr::filter(!is.na(values)) |>
      dplyr::summarise(values = sum(values)) |>
      dplyr::ungroup()

    final_table <-
      dplyr::bind_rows(
        final_table_1,
        final_table_2,
        final_table_3,
        final_table_4
      ) |>
      # dplyr::filter(!grepl(pattern = "^Other", x = ids)) |>
      dplyr::filter(!is.na(sample))

    if (type == "analysis") {
      if (rescale == TRUE) {
        final_table <- final_table |>
          dplyr::group_by(sample, species) |>
          dplyr::mutate(values = values / (sum(values) / 3)) |> #' because 3 levels
          dplyr::ungroup()
      } else {
      final_table <- final_table |>
        dplyr::group_by(sample) |>
        dplyr::mutate(values = values / (sum(values) / 3)) |> #' because 3 levels
        dplyr::ungroup()
      }
    }

    return(final_table)
  }
