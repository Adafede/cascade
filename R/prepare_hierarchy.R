#' Prepare hierarchy
#'
#' @param dataframe Dataframe
#' @param type Type
#' @param detector Detector
#' @param rescale Rescale
#'
#' @return A dataframe with prepared hierarchy
#'
#' @examples NULL
prepare_hierarchy <-
  function(dataframe, type = "analysis", detector = "ms", rescale = FALSE) {
    stopifnot(
      "'type' must be either 'analysis' or 'literature'" = type %in%
        c("analysis", "literature")
    )

    if (type == "analysis") {
      df <- dataframe |>
        tidytable::mutate(
          best_candidate_1 = tidytable::if_else(
            condition = !is.na(smiles_2D),
            true = best_candidate_1,
            false = "Other"
          )
        ) |>
        tidytable::mutate(
          best_candidate_2 = tidytable::if_else(
            condition = !is.na(smiles_2D),
            true = best_candidate_2,
            false = paste(best_candidate_1, "notAnnotated", sep = " ")
          )
        ) |>
        tidytable::mutate(
          best_candidate_3 = tidytable::if_else(
            condition = !is.na(smiles_2D),
            true = best_candidate_3,
            false = paste(best_candidate_2, "notAnnotated", sep = " ")
          )
        ) |>
        tidytable::mutate(
          best_candidate_2 = tidytable::if_else(
            condition = best_candidate_2 != "notClassified",
            true = best_candidate_2,
            false = paste(best_candidate_1, "notClassified", sep = " ")
          )
        ) |>
        tidytable::mutate(
          best_candidate_3 = tidytable::if_else(
            condition = best_candidate_3 != "notClassified",
            true = best_candidate_3,
            false = paste(best_candidate_2, "notClassified", sep = " ")
          )
        ) |>
        tidytable::mutate(
          best_candidate_1 = tidytable::if_else(
            condition = best_candidate_3 !=
              "notClassified notClassified notClassified",
            true = best_candidate_1,
            false = "Other"
          )
        ) |>
        tidytable::mutate(
          best_candidate_2 = tidytable::if_else(
            condition = best_candidate_3 !=
              "notClassified notClassified notClassified",
            true = best_candidate_2,
            false = "Other notClassified"
          )
        ) |>
        tidytable::mutate(
          best_candidate_3 = tidytable::if_else(
            condition = best_candidate_3 !=
              "notClassified notClassified notClassified",
            true = best_candidate_3,
            false = "Other notClassified notClassified"
          )
        ) |>
        tidytable::mutate(chemical_pathway = best_candidate_1)

      df_notConfident <- df |>
        tidytable::filter(grepl(
          pattern = "notConfident",
          x = best_candidate_2
        )) |>
        tidytable::mutate(
          smiles_2D = NA,
          best_candidate_1 = "Other",
          best_candidate_2 = "Other notConfident",
          best_candidate_3 = "Other notConfident notConfident",
          inchikey_2D = NA
        ) |>
        tidytable::arrange(tidytable::desc(intensity)) |>
        tidytable::arrange(tidytable::desc(comparison_score)) |>
        tidytable::distinct(id, peak_id, sample, .keep_all = TRUE)

      df_confident <- df |>
        tidytable::filter(
          !grepl(pattern = "notConfident", x = best_candidate_2)
        )

      dataframe2 <- tidytable::bind_rows(df_confident, df_notConfident) |>
        tidytable::group_by(chemical_pathway)
    } else {
      dataframe2 <- dataframe |>
        tidytable::mutate(
          best_candidate_1 = tidytable::if_else(
            condition = !is.na(best_candidate_1),
            true = best_candidate_1,
            false = "Other"
          ),
          best_candidate_2 = tidytable::if_else(
            condition = !is.na(best_candidate_2),
            true = best_candidate_2,
            false = paste(best_candidate_1, "notClassified", sep = " ")
          ),
          best_candidate_3 = tidytable::if_else(
            condition = !is.na(best_candidate_3),
            true = best_candidate_3,
            false = paste(best_candidate_2, "notClassified", sep = " ")
          )
        ) |>
        tidytable::mutate(
          chemical_pathway = tidytable::if_else(
            condition = grepl(pattern = "not", x = best_candidate_1),
            true = "Other",
            false = best_candidate_1
          )
        )
    }

    ## dirty temp fix because of bad npclassifier naming
    dataframe2 <- dataframe2 |>
      tidytable::rowwise() |>
      tidytable::mutate(
        safety = tidytable::if_else(
          condition = best_candidate_3 == best_candidate_2,
          true = "Y",
          false = "N"
        )
      ) |>
      tidytable::mutate(
        best_candidate_3 = tidytable::if_else(
          condition = safety == "Y",
          true = paste(best_candidate_3, "sub"),
          false = best_candidate_3
        )
      )

    parents <- dataframe2 |>
      tidytable::distinct(
        labels = best_candidate_1,
        ids = best_candidate_1,
        chemical_pathway
      ) |>
      tidytable::mutate(parents = "") |>
      tidytable::distinct() |>
      tidytable::mutate(
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
      tidytable::mutate(
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
      tidytable::mutate(
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
      tidytable::distinct()

    children_1 <- dataframe2 |>
      tidytable::distinct(
        best_candidate_1,
        best_candidate_2,
        chemical_pathway
      ) |>
      tidytable::mutate(
        ids = paste(best_candidate_1, best_candidate_2, sep = "-")
      ) |>
      tidytable::distinct(labels = best_candidate_2, ids, chemical_pathway) |>
      tidytable::mutate(
        parents = gsub(
          pattern = "-.*",
          replacement = "",
          x = ids
        )
      ) |>
      tidytable::mutate(tidytable::across(
        .cols = tidytable::everything(),
        .fns = ~ gsub(
          x = .x,
          pattern = "-NA$",
          replacement = "",
        )
      )) |>
      tidytable::filter(!is.na(labels)) |>
      tidytable::distinct()

    children_2 <- dataframe2 |>
      tidytable::distinct(
        best_candidate_2,
        best_candidate_3,
        chemical_pathway
      ) |>
      tidytable::mutate(
        ids = paste(best_candidate_2, best_candidate_3, sep = "-")
      ) |>
      tidytable::distinct(
        labels = best_candidate_3,
        ids,
        join = best_candidate_2,
        chemical_pathway
      ) |>
      tidytable::mutate(
        join = gsub(
          pattern = "([a-zA-Z]|\\))(-)([a-zA-Z]|[0-9])(.*)",
          replacement = "\\1",
          x = join
        )
      ) |>
      tidytable::mutate(tidytable::across(
        .cols = tidytable::everything(),
        .fns = ~ gsub(
          x = .x,
          pattern = "-NA$",
          replacement = "",
        )
      )) |>
      tidytable::mutate(tidytable::across(
        .cols = tidytable::everything(),
        .fns = ~ gsub(
          x = .x,
          pattern = "^NA$",
          replacement = NA,
        )
      )) |>
      tidytable::full_join(
        children_1,
        by = c("join" = "labels", "chemical_pathway" = "chemical_pathway")
      ) |>
      tidytable::distinct(
        ids = ids.x,
        labels,
        parents = ids.y,
        chemical_pathway
      ) |>
      tidytable::filter(!is.na(labels)) |>
      tidytable::distinct()

    genealogy <- tidytable::bind_rows(parents, children_1, children_2) |>
      tidytable::ungroup() |>
      tidytable::distinct() |>
      tidytable::group_by(ids) |>
      tidytable::add_count() |> ## for ambiguous classes
      tidytable::filter(
        parents != "" |
          !is.na(ids) |
          !is.na(labels)
      ) |>
      tidytable::ungroup()

    suppressMessages(
      table <- dataframe2 |>
        tidytable::mutate(
          ids = paste0(
            best_candidate_2,
            "-",
            best_candidate_3
          )
        ) |>
        tidytable::mutate(
          ids = gsub(
            x = ids,
            pattern = "-NA$",
            replacement = "",
          )
        ) |>
        tidytable::mutate(
          ids = gsub(
            x = ids,
            pattern = "^NA$",
            replacement = NA,
          )
        ) |>
        tidytable::full_join(genealogy) |>
        tidytable::mutate(
          best_candidate_3 = tidytable::if_else(
            condition = is.na(best_candidate_3) & !is.na(best_candidate_2),
            false = best_candidate_3,
            true = paste0(best_candidate_1, "-", best_candidate_2)
          )
        )
    )

    if (type == "analysis") {
      table <- table |>
        tidytable::distinct(
          feature_id,
          smiles_2D,
          inchikey_2D,
          score_biological,
          score_chemical,
          score_final,
          best_candidate_organism,
          labels,
          ids,
          chemical_pathway,
          parents,
          sample,
          intensity,
          species,
          n
        )

      table_1 <- table |>
        tidytable::group_by(parents, ids, labels, sample) |>
        tidytable::add_count(name = "values") |>
        tidytable::select(
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
        tidytable::distinct()

      table_1_1 <- table_1 |>
        tidytable::group_by(chemical_pathway, parents, ids, sample) |>
        tidytable::mutate(values_2 = sum(values)) |>
        tidytable::ungroup()
    } else {
      table_1 <- table |>
        tidytable::group_by(parents, ids, labels, organism) |>
        tidytable::add_count(name = "values") |>
        tidytable::ungroup() |>
        tidytable::select(
          chemical_pathway,
          parents,
          ids,
          labels,
          values,
          organism,
          sample,
          n
        ) |>
        tidytable::distinct()

      table_1_1 <- table_1 |>
        tidytable::group_by(chemical_pathway, parents, ids, organism) |>
        tidytable::mutate(values_2 = sum(values)) |>
        tidytable::ungroup()
    }

    top_parents_table <- table_1_1 |>
      tidytable::distinct(
        chemical_pathway,
        parents,
        ids,
        labels,
        values_2,
        n
      ) |>
      tidytable::inner_join(
        table_1_1 |>
          tidytable::distinct(
            chemical_pathway,
            parents,
            ids,
            labels,
            values_2,
            n
          ),
        by = c("ids" = "parents")
      ) |>
      tidytable::filter(!grepl(pattern = "-", x = parents)) |>
      tidytable::distinct(
        chemical_pathway = chemical_pathway.x,
        parents,
        ids,
        labels = labels.x,
        values_2 = values_2.y,
        n.x
      ) |>
      tidytable::group_by(chemical_pathway, parents) |>
      tidytable::mutate(values_3 = sum(!is.na(values_2)) / as.numeric(n.x)) |>
      tidytable::filter(parents != "") |>
      tidytable::distinct(chemical_pathway, parents, values_3) |>
      tidytable::ungroup() |>
      tidytable::top_n(n = length(microshades), wt = values_3) |>
      tidytable::arrange(tidytable::desc(values_3)) |>
      tidytable::select(-values_3) |>
      tidytable::mutate(ids = parents, labels = parents) |>
      tidytable::mutate(parents = "", new_labels = labels) |>
      tidytable::distinct()

    if (nrow(top_parents_table > length(microshades))) {
      top_parents_table <- top_parents_table |>
        utils::head(length(microshades)) ## in case of equal numbers among classes
    }

    top_parents <- top_parents_table$labels

    table_2 <- table_1_1 |>
      tidytable::group_by(chemical_pathway, parents) |>
      tidytable::summarize(sum = sum(values), .groups = "drop") |>
      tidytable::ungroup()

    suppressMessages(
      table_3 <-
        tidytable::inner_join(table_2, children_1, by = c("parents" = "ids")) |>
        tidytable::rename(chemical_pathway = chemical_pathway.y) |>
        tidytable::group_by(chemical_pathway, i.parents) |>
        tidytable::mutate(n = sum(sum)) |>
        tidytable::ungroup() |>
        tidytable::mutate(
          n = tidytable::if_else(
            condition = i.parents == "Other",
            true = n + 666,
            false = n + 0
          )
        ) |>
        tidytable::select(
          chemical_pathway,
          parents = i.parents,
          ids = parents,
          labels,
          sum,
          n
        )
    )

    table_3_1 <- table_3 |>
      tidytable::distinct(n, .keep_all = TRUE) |>
      ## COMMENT replacement because not working
      tidytable::arrange(tidytable::desc(n)) |>
      tidytable::mutate(row_number = tidytable::row_number()) |>
      tidytable::filter(row_number <= 4L) |>
      tidytable::select(-row_number) |>
      ## COMMENT replacement because not working
      # tidytable::slice_max(order_by = n, n = 4L, with_ties = FALSE) |>
      tidytable::select(chemical_pathway, parents, ids, labels)

    top_medium_table <- table_3 |>
      tidytable::filter(parents %in% top_parents) |>
      tidytable::group_by(parents) |>
      tidytable::distinct(sum, .keep_all = TRUE) |>
      ## COMMENT replacement because not working
      tidytable::arrange(tidytable::desc(sum)) |>
      tidytable::mutate(row_number = tidytable::row_number(), .by = parents) |>
      tidytable::filter(row_number <= 4L) |>
      tidytable::select(-row_number) |>
      ## COMMENT replacement because not working
      # tidytable::slice_max(order_by = sum, by = parents, n = 4L, with_ties = FALSE) |>
      tidytable::select(chemical_pathway, parents, ids, labels) |>
      tidytable::mutate(new_labels = labels) |>
      tidytable::ungroup() |>
      tidytable::distinct()

    # top_medium <- top_medium_table$ids

    suppressMessages(
      low_medium_table <-
        tidytable::anti_join(
          table_3,
          tidytable::bind_rows(table_3_1, top_medium_table)
        ) |>
        tidytable::select(chemical_pathway, parents, ids, labels) |>
        tidytable::filter(!is.na(parents)) |>
        tidytable::distinct()
    )

    low_medium_table_1_a <- low_medium_table |>
      tidytable::filter(parents %in% top_parents) |>
      tidytable::mutate(new_labels = paste(parents, "Other")) |>
      tidytable::mutate(ids = paste(parents, new_labels, sep = "-")) |>
      tidytable::mutate(labels = new_labels) |>
      tidytable::distinct()

    low_medium_table_1_b <- low_medium_table |>
      tidytable::filter(parents %in% top_parents) |>
      tidytable::mutate(new_labels = paste(parents, "Other")) |>
      tidytable::mutate(
        parents = paste(parents, new_labels, sep = "-"),
        ids = paste(new_labels, labels, sep = "-")
      ) |>
      tidytable::mutate(new_labels = labels) |>
      tidytable::distinct()

    low_medium_table_1 <-
      tidytable::bind_rows(low_medium_table_1_a, low_medium_table_1_b) |>
      tidytable::distinct()

    low_medium_table_2 <- low_medium_table |>
      tidytable::filter(!parents %in% top_parents) |>
      tidytable::mutate(new_labels = paste("Other", parents)) |>
      tidytable::mutate(parents = "Other") |>
      tidytable::mutate(ids = paste(parents, new_labels, sep = "-")) |>
      tidytable::distinct()

    genealogy_new_med <-
      tidytable::bind_rows(
        top_parents_table,
        top_medium_table,
        low_medium_table_2,
        low_medium_table_1
      ) |>
      tidytable::distinct()

    table_next <- genealogy |>
      tidytable::filter(!labels %in% genealogy_new_med$labels) |>
      tidytable::filter(parents %in% top_medium_table$ids) |>
      tidytable::filter(parents != "") |>
      tidytable::mutate(
        new_labels = paste(
          gsub(
            pattern = "^([^-]*)-(.*)",
            replacement = "",
            x = parents
          ),
          labels
        ) |>
          trimws()
      ) |>
      tidytable::mutate(
        ids = paste(
          gsub(
            pattern = "^([^-]*)-(.*)",
            replacement = "\\2",
            x = parents
          ),
          new_labels,
          sep = "-"
        )
      ) |>
      tidytable::mutate(
        new_labels = tidytable::if_else(
          condition = grepl(pattern = "^Other", x = parents),
          true = labels,
          false = new_labels
        )
      )

    genealogy_new_med_2 <-
      tidytable::bind_rows(genealogy_new_med, table_next) |>
      tidytable::ungroup() |>
      tidytable::distinct()

    suppressMessages(
      table_children <- genealogy |>
        tidytable::filter(!labels %in% genealogy_new_med$ids) |>
        tidytable::select(labels) |>
        tidytable::left_join(genealogy_new_med_2) |>
        tidytable::filter(grepl(pattern = "-", x = parents)) |>
        tidytable::inner_join(genealogy, by = c("parents" = "ids")) |>
        tidytable::select(
          chemical_pathway = chemical_pathway.x,
          parents,
          ids,
          labels = labels.x,
          new_labels
        )
    )

    table_other <- genealogy_new_med_2 |>
      tidytable::filter(grepl(pattern = " Other-", x = ids)) |>
      tidytable::distinct(chemical_pathway, ids = parents) |>
      tidytable::mutate(
        parents = (gsub(
          pattern = "-.*",
          replacement = "\\1",
          x = ids
        ))
      ) |>
      tidytable::mutate(
        labels = (gsub(
          pattern = ".*-",
          replacement = "\\1",
          x = ids
        )),
        new_labels = labels
      )

    genealogy_new_med_3 <-
      tidytable::bind_rows(genealogy_new_med_2, table_children, table_other) |>
      tidytable::distinct()

    suppressMessages(
      missing_children_1 <- genealogy_new_med_3 |>
        tidytable::filter(grepl(pattern = " Other$", x = parents)) |>
        tidytable::distinct(
          chemical_pathway,
          parents = ids,
          best_candidate_2 = labels
        ) |>
        tidytable::left_join(
          dataframe2 |>
            tidytable::mutate(
              best_candidate_3 = as.character(best_candidate_3),
              best_candidate_2 = as.character(best_candidate_2)
            ) |>
            tidytable::distinct(best_candidate_3, best_candidate_2)
        ) |>
        tidytable::mutate(
          ids = paste(best_candidate_2, best_candidate_3, sep = "-"),
          labels = best_candidate_3,
          new_labels = best_candidate_3
        ) |>
        tidytable::select(chemical_pathway, parents, ids, labels, new_labels) |>
        tidytable::mutate(tidytable::across(
          .cols = tidytable::everything(),
          .fns = ~ gsub(
            x = .x,
            pattern = "-NA$",
            replacement = "",
          )
        ))
    )

    suppressMessages(
      missing_children_2 <- genealogy_new_med_3 |>
        tidytable::filter(grepl(pattern = "^Other", x = parents)) |>
        tidytable::distinct(
          chemical_pathway,
          parents = ids,
          best_candidate_2 = labels
        ) |>
        tidytable::left_join(
          dataframe2 |>
            tidytable::mutate(
              best_candidate_3 = as.character(best_candidate_3),
              best_candidate_2 = as.character(best_candidate_2)
            ) |>
            tidytable::distinct(best_candidate_3, best_candidate_2)
        ) |>
        tidytable::mutate(
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
        tidytable::select(chemical_pathway, parents, ids, labels, new_labels) |>
        tidytable::distinct() |>
        tidytable::mutate(tidytable::across(
          .cols = tidytable::everything(),
          .fns = ~ gsub(
            x = .x,
            pattern = "-NA$",
            replacement = "",
          )
        )) |>
        tidytable::filter(!is.na(labels))
    )

    missing_children <-
      tidytable::bind_rows(missing_children_1, missing_children_2)

    genealogy_new_med_4 <-
      tidytable::bind_rows(
        genealogy_new_med_3 |>
          tidytable::mutate(chemical_pathway = as.character(chemical_pathway)),
        missing_children
      ) |>
      tidytable::distinct(chemical_pathway, parents, ids, labels, new_labels) |>
      tidytable::mutate(link = ids |> gsub(pattern = "-.*", replacement = ""))

    table_new <- dataframe2 |>
      tidytable::filter(!is.na(species)) |>
      tidytable::mutate(
        best_candidate_3 = as.character(best_candidate_3),
        chemical_pathway = as.character(chemical_pathway)
      ) |>
      tidytable::full_join(
        genealogy_new_med_4,
        by = c(
          "best_candidate_2" = "link",
          "best_candidate_3" = "labels",
          "chemical_pathway" = "chemical_pathway"
        )
      )

    if (type == "analysis") {
      if (nrow(table_new) > 0) {
        table_new <- table_new |>
          tidytable::rowwise() |>
          tidytable::filter(
            id %in%
              sample |
              ## Comment these if
              is.na(id)
          ) |> ## using wrong feature table
          tidytable::ungroup()
      } else {
        table_new <- table_new |>
          tidytable::mutate(integral = NA, new_labels = NA)
      }
      table_new <- table_new |>
        tidytable::mutate(
          intensity = switch(detector, "ms" = intensity, "cad" = integral)
        ) |>
        tidytable::distinct(
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
        tidytable::group_by(parents, ids, sample, intensity) |>
        tidytable::add_count(name = "values") |>
        tidytable::ungroup() |>
        tidytable::select(
          parents,
          ids,
          labels,
          values,
          sample,
          intensity,
          species
        ) |>
        tidytable::distinct() |>
        tidytable::filter(!is.na(ids))
    } else {
      table_1_new <- table_new |>
        tidytable::mutate(labels = best_candidate_3) |>
        tidytable::group_by(parents, ids, sample) |>
        tidytable::add_count(name = "values") |>
        tidytable::ungroup() |>
        tidytable::select(
          parents,
          ids,
          labels,
          values,
          sample,
          organism,
          species
        ) |>
        tidytable::distinct() |>
        tidytable::filter(!is.na(ids))
    }

    additional_row <- table_1_new |>
      tidytable::filter(labels == "Other other") |>
      tidytable::mutate(
        parents = "Other-Other other",
        ids = "Other other - Other other other",
        labels = "Other other other"
      )

    table_1_1_new <- tidytable::bind_rows(
      table_1_new |>
        tidytable::mutate(
          labels = labels |>
            as.character()
        ),
      additional_row
    ) |>
      tidytable::ungroup()

    final_table_3_1 <- table_1_1_new |>
      tidytable::filter(ids %in% low_medium_table_1_a$ids) |>
      tidytable::distinct(parents, ids, labels)

    final_table_4_1 <- table_1_1_new |>
      tidytable::filter(
        grepl(pattern = "-", x = parents) &
          parents %in% missing_children_1$parents
      ) |>
      tidytable::group_by(
        c("parents", "ids", "labels", "sample", "species")
      ) |>
      tidytable::summarize(
        values = sum(switch(
          type,
          "analysis" = intensity,
          "literature" = values
        )),
        .groups = "drop"
      ) |>
      tidytable::filter(!is.na(labels)) |>
      tidytable::mutate(
        join = gsub(
          pattern = "-.*",
          replacement = "",
          x = parents
        )
      ) |>
      tidytable::ungroup()

    final_table_4 <- final_table_4_1 |>
      tidytable::select(-join)

    final_table_3_2 <- final_table_3_1 |>
      tidytable::inner_join(final_table_4_1, by = c("labels" = "join")) |>
      tidytable::select(
        parents = ids.x,
        ids = parents.y,
        labels = ids.y,
        sample,
        species,
        values
      ) |>
      tidytable::mutate(
        labels = gsub(
          pattern = "([a-zA-Z]|\\))(-)([a-zA-Z]|[0-9])(.*)",
          replacement = "\\1",
          x = labels
        )
      ) |>
      tidytable::group_by(
        c("parents", "ids", "labels", "sample", "species")
      ) |>
      tidytable::summarize(values = sum(values), .groups = "drop") |>
      tidytable::ungroup()

    final_table_3_3 <- table_1_1_new |>
      tidytable::filter(!is.na(species)) |>
      tidytable::filter(
        grepl(pattern = "-", x = parents) &
          !parents %in% missing_children_1$parents
      ) |>
      tidytable::group_by(c("parents", "ids", "labels", "sample", "species")) |>
      tidytable::summarize(
        values = sum(switch(
          type,
          "analysis" = intensity,
          "literature" = values
        )),
        .groups = "drop"
      ) |>
      tidytable::ungroup()

    final_table_3 <-
      tidytable::bind_rows(final_table_3_3, final_table_3_2) |>
      tidytable::distinct() |>
      tidytable::arrange(parents, ids)

    final_table_2 <- table_1_1_new |>
      tidytable::select(
        -tidytable::any_of(c(
          "values",
          "sample",
          "intensity",
          "species",
          "organism"
        ))
      ) |>
      tidytable::left_join(
        final_table_3 |>
          tidytable::select(tidytable::any_of(
            c("values", "sample", "parents", "species")
          )),
        by = c("ids" = "parents")
      ) |>
      tidytable::group_by(
        c("parents", "ids", "labels", "sample", "species")
      ) |>
      tidytable::summarize(values = sum(values), .groups = "drop") |>
      tidytable::ungroup() |>
      tidytable::filter(!is.na(values))

    final_table_1 <- table_1_1_new |>
      tidytable::select(
        -tidytable::any_of(c(
          "values",
          "sample",
          "intensity",
          "species",
          "organism"
        ))
      ) |>
      tidytable::left_join(
        final_table_2 |>
          tidytable::select(tidytable::any_of(
            c("values", "sample", "parents", "species")
          )),
        by = c("ids" = "parents")
      ) |>
      tidytable::group_by(
        c("parents", "ids", "labels", "sample", "species")
      ) |>
      tidytable::filter(!is.na(values)) |>
      tidytable::summarize(values = sum(values), .groups = "drop") |>
      tidytable::ungroup()

    final_table <-
      tidytable::bind_rows(
        final_table_1,
        final_table_2,
        final_table_3,
        final_table_4
      ) |>
      # tidytable::filter(!grepl(pattern = "^Other", x = ids)) |>
      tidytable::filter(!is.na(sample))

    if (type == "analysis") {
      if (rescale == TRUE) {
        final_table <- final_table |>
          tidytable::group_by(sample, species) |>
          tidytable::mutate(values = values / (sum(values) / 3)) |> ## because 3 levels
          tidytable::ungroup()
      } else {
        final_table <- final_table |>
          tidytable::group_by(sample) |>
          tidytable::mutate(values = values / 3) |> ## because 3 levels
          tidytable::ungroup()
      }
    }

    return(final_table)
  }
