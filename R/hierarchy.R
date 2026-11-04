#' Reorder factor levels by two variables (base R replacement for forcats::fct_reorder2)
#'
#' @param f Factor or character vector.
#' @param x First sorting vector.
#' @param y Second sorting vector.
#' @param desc Logical; sort in descending order (default `TRUE`).
#'
#' @return A factor with reordered levels.
#' @keywords internal
#' @noRd
fct_reorder2_base <- function(f, x, y, desc = TRUE) {
  f <- as.factor(f)
  idx <- split(seq_along(f), f, drop = FALSE)
  last_i <- vapply(
    idx,
    function(i) {
      if (length(i) == 0L) {
        return(NA_integer_)
      }
      ord <- order(x[i], y[i], na.last = FALSE)
      i[ord[length(ord)]]
    },
    integer(1L)
  )
  keys <- y[last_i]
  factor(f, levels = levels(f)[order(keys, decreasing = desc, na.last = TRUE)])
}

#' Dense rank of a vector (base R replacement for tidytable::dense_rank)
#'
#' @param x Vector to rank.
#'
#' @return Integer vector of dense ranks (`1..K`).
#' @keywords internal
#' @noRd
dense_rank_base <- function(x) {
  match(x, sort(unique(x)))
}

#' Prepare chemical class hierarchy for sunburst/treemap and histogram plots
#'
#' @param dataframe Data frame or `data.table` of features/compounds.
#' @param type Either `"analysis"` or `"literature"`.
#' @param detector Detector type (`"ms"` or `"cad"`).
#' @param rescale Logical; whether to rescale values to 100% per sample.
#'
#' @return A `data.table` with hierarchical `parents`, `ids`, `labels`, and
#'   `values`.
#' @keywords internal
#' @noRd
prepare_hierarchy <- function(
  dataframe,
  type = "analysis",
  detector = "ms",
  rescale = FALSE
) {
  stopifnot(
    "'type' must be either 'analysis' or 'literature'" = type %in%
      c("analysis", "literature")
  )

  dt_in <- data.table::copy(data.table::as.data.table(dataframe))

  if (type == "analysis") {
    has_smi <- !is.na(dt_in$smiles_2D)
    dt_in[, best_candidate_1 := ifelse(has_smi, best_candidate_1, "Other")]
    dt_in[,
      best_candidate_2 := ifelse(
        has_smi,
        best_candidate_2,
        paste(best_candidate_1, "notAnnotated", sep = " ")
      )
    ]
    dt_in[,
      best_candidate_3 := ifelse(
        has_smi,
        best_candidate_3,
        paste(best_candidate_2, "notAnnotated", sep = " ")
      )
    ]
    dt_in[,
      best_candidate_2 := ifelse(
        !is.na(best_candidate_2) & best_candidate_2 != "notClassified",
        best_candidate_2,
        paste(best_candidate_1, "notClassified", sep = " ")
      )
    ]
    dt_in[,
      best_candidate_3 := ifelse(
        !is.na(best_candidate_3) & best_candidate_3 != "notClassified",
        best_candidate_3,
        paste(best_candidate_2, "notClassified", sep = " ")
      )
    ]
    is_unclass3 <- !is.na(dt_in$best_candidate_3) &
      dt_in$best_candidate_3 == "notClassified notClassified notClassified"
    dt_in[
      is_unclass3,
      `:=`(
        best_candidate_1 = "Other",
        best_candidate_2 = "Other notClassified",
        best_candidate_3 = "Other notClassified notClassified"
      )
    ]
    dt_in[, chemical_pathway := best_candidate_1]

    is_not_conf <- grepl("notConfident", dt_in$best_candidate_2)
    df_not_conf <- data.table::copy(dt_in[is_not_conf])
    if (nrow(df_not_conf) > 0L) {
      df_not_conf[, `:=`(
        smiles_2D = NA_character_,
        best_candidate_1 = "Other",
        best_candidate_2 = "Other notConfident",
        best_candidate_3 = "Other notConfident notConfident",
        inchikey_2D = NA_character_
      )]
      data.table::setorder(
        df_not_conf,
        -comparison_score,
        -intensity,
        na.last = TRUE
      )
      df_not_conf <- unique(df_not_conf, by = c("id", "peak_id", "sample"))
    }
    df_conf <- dt_in[!is_not_conf]
    dataframe2 <- data.table::rbindlist(
      list(df_conf, df_not_conf),
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    dataframe2 <- dt_in
    dataframe2[is.na(best_candidate_1), best_candidate_1 := "Other"]
    dataframe2[
      is.na(best_candidate_2),
      best_candidate_2 := paste(best_candidate_1, "notClassified", sep = " ")
    ]
    dataframe2[
      is.na(best_candidate_3),
      best_candidate_3 := paste(best_candidate_2, "notClassified", sep = " ")
    ]
    dataframe2[,
      chemical_pathway := ifelse(
        grepl("not", best_candidate_1),
        "Other",
        best_candidate_1
      )
    ]
  }

  same_23 <- !is.na(dataframe2$best_candidate_3) &
    !is.na(dataframe2$best_candidate_2) &
    dataframe2$best_candidate_3 == dataframe2$best_candidate_2
  dataframe2[, safety := ifelse(same_23, "Y", "N")]
  dataframe2[safety == "Y", best_candidate_3 := paste(best_candidate_3, "sub")]

  clean_other <- function(x) {
    x <- gsub("notAnnotated", "Other", x, fixed = TRUE)
    x <- gsub("notClassified", "Other", x, fixed = TRUE)
    gsub("notConfident", "Other", x, fixed = TRUE)
  }

  parents <- unique(dataframe2[, list(
    labels = clean_other(best_candidate_1),
    ids = clean_other(best_candidate_1),
    chemical_pathway = chemical_pathway,
    parents = ""
  )])

  children_1 <- unique(dataframe2[, list(
    best_candidate_1,
    best_candidate_2,
    chemical_pathway
  )])
  children_1[, ids := paste(best_candidate_1, best_candidate_2, sep = "-")]
  children_1 <- unique(children_1[, list(
    labels = best_candidate_2,
    ids,
    chemical_pathway,
    parents = gsub("-.*", "", ids)
  )])
  for (col in names(children_1)) {
    data.table::set(
      children_1,
      j = col,
      value = gsub("-NA$", "", children_1[[col]])
    )
  }
  children_1 <- unique(children_1[!is.na(labels)])

  children_2_pre <- unique(dataframe2[, list(
    best_candidate_2,
    best_candidate_3,
    chemical_pathway
  )])
  children_2_pre[, ids := paste(best_candidate_2, best_candidate_3, sep = "-")]
  children_2_pre <- unique(children_2_pre[, list(
    labels = best_candidate_3,
    ids,
    join = gsub(
      "([a-zA-Z]|\\))(-)([a-zA-Z]|[0-9])(.*)",
      "\\1",
      best_candidate_2
    ),
    chemical_pathway
  )])
  for (col in names(children_2_pre)) {
    val <- gsub("-NA$", "", children_2_pre[[col]])
    val[val == "NA"] <- NA_character_
    data.table::set(children_2_pre, j = col, value = val)
  }

  c2_joined <- dt_full_join(
    children_2_pre,
    children_1,
    by = c("join" = "labels", "chemical_pathway" = "chemical_pathway")
  )
  children_2 <- unique(c2_joined[
    !is.na(labels),
    list(
      ids = ids.x,
      labels = labels,
      parents = ids.y,
      chemical_pathway = chemical_pathway
    )
  ])

  genealogy <- unique(data.table::rbindlist(
    list(parents, children_1, children_2),
    use.names = TRUE,
    fill = TRUE
  ))
  genealogy[, n := .N, by = "ids"]
  genealogy <- genealogy[parents != "" | !is.na(ids) | !is.na(labels)]

  df2_ids <- data.table::copy(dataframe2)
  df2_ids[,
    ids := gsub("-NA$", "", paste0(best_candidate_2, "-", best_candidate_3))
  ]
  df2_ids[ids == "NA", ids := NA_character_]

  table_joined <- dt_full_join(df2_ids, genealogy)
  table_joined[
    is.na(best_candidate_3) & !is.na(best_candidate_2),
    best_candidate_3 := paste0(best_candidate_1, "-", best_candidate_2)
  ]

  if (type == "analysis") {
    table_sub <- unique(table_joined[, list(
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
    )])
    table_sub[, values := .N, by = c("parents", "ids", "labels", "sample")]
    table_1 <- unique(table_sub[, list(
      chemical_pathway,
      parents,
      ids,
      labels,
      values,
      sample,
      intensity,
      species,
      n
    )])
    table_1_1 <- data.table::copy(table_1)
    table_1_1[,
      values_2 := sum(values),
      by = c("chemical_pathway", "parents", "ids", "sample")
    ]
  } else {
    table_joined[, values := .N, by = c("parents", "ids", "labels", "organism")]
    table_1 <- unique(table_joined[, list(
      chemical_pathway,
      parents,
      ids,
      labels,
      values,
      organism,
      sample,
      n
    )])
    table_1_1 <- data.table::copy(table_1)
    table_1_1[,
      values_2 := sum(values),
      by = c("chemical_pathway", "parents", "ids", "organism")
    ]
  }

  t1_dist <- unique(table_1_1[, list(
    chemical_pathway,
    parents,
    ids,
    labels,
    values_2,
    n
  )])
  tp_join <- dt_inner_join(
    t1_dist,
    t1_dist,
    by = c("ids" = "parents")
  )
  tp_join <- unique(tp_join[
    !grepl("-", parents),
    list(
      chemical_pathway = chemical_pathway.x,
      parents = parents,
      ids = ids,
      labels = labels.x,
      values_2 = values_2.y,
      n.x = n.x
    )
  ])
  tp_join[,
    values_3 := sum(!is.na(values_2)) / as.numeric(n.x),
    by = c("chemical_pathway", "parents")
  ]
  tp_join <- unique(tp_join[
    parents != "",
    list(chemical_pathway, parents, values_3)
  ])
  data.table::setorder(tp_join, -values_3, na.last = TRUE)
  if (nrow(tp_join) > length(microshades)) {
    tp_join <- utils::head(tp_join, length(microshades))
  }
  top_parents_table <- unique(tp_join[, list(
    chemical_pathway = chemical_pathway,
    parents = "",
    ids = parents,
    labels = parents,
    new_labels = parents
  )])

  top_parents <- top_parents_table$labels

  table_2 <- table_1_1[,
    list(sum = sum(values)),
    keyby = c("chemical_pathway", "parents")
  ]
  data.table::setkey(table_2, NULL)

  t3_join <- dt_inner_join(
    table_2,
    children_1,
    by = c("parents" = "ids")
  )
  t3_join[, n := sum(sum), by = c("chemical_pathway.y", "i.parents")]
  t3_join[i.parents == "Other", n := n + 666]
  table_3 <- t3_join[, list(
    chemical_pathway = chemical_pathway.y,
    parents = i.parents,
    ids = parents,
    labels = labels,
    sum = sum,
    n = n
  )]

  t3_1 <- unique(table_3, by = "n")
  data.table::setorder(t3_1, -n, na.last = TRUE)
  table_3_1 <- utils::head(t3_1, 4L)[, list(
    chemical_pathway,
    parents,
    ids,
    labels
  )]

  tm <- table_3[parents %in% top_parents]
  tm <- unique(tm, by = c("parents", "sum"))
  data.table::setorder(tm, -sum, na.last = TRUE)
  tm[, row_number := seq_len(.N), by = "parents"]
  top_medium_table <- unique(tm[
    row_number <= 4L,
    list(chemical_pathway, parents, ids, labels, new_labels = labels)
  ])

  excl <- data.table::rbindlist(
    list(
      table_3_1,
      top_medium_table[, list(chemical_pathway, parents, ids, labels)]
    ),
    use.names = TRUE,
    fill = TRUE
  )
  low_medium_table <- unique(
    table_3[!excl, on = c("chemical_pathway", "parents", "ids", "labels")][
      !is.na(parents),
      list(chemical_pathway, parents, ids, labels)
    ]
  )

  lm_top <- low_medium_table[parents %in% top_parents]
  low_medium_table_1_a <- unique(lm_top[, list(
    chemical_pathway = chemical_pathway,
    parents = parents,
    ids = paste(parents, paste(parents, "Other"), sep = "-"),
    labels = paste(parents, "Other"),
    new_labels = paste(parents, "Other")
  )])
  low_medium_table_1_b <- unique(lm_top[, list(
    chemical_pathway = chemical_pathway,
    parents = paste(parents, paste(parents, "Other"), sep = "-"),
    ids = paste(paste(parents, "Other"), labels, sep = "-"),
    labels = labels,
    new_labels = labels
  )])
  low_medium_table_1 <- unique(data.table::rbindlist(
    list(low_medium_table_1_a, low_medium_table_1_b),
    use.names = TRUE
  ))

  lm_not_top <- low_medium_table[!parents %in% top_parents]
  low_medium_table_2 <- unique(lm_not_top[, list(
    chemical_pathway = chemical_pathway,
    parents = "Other",
    ids = paste("Other", paste("Other", parents), sep = "-"),
    labels = labels,
    new_labels = paste("Other", parents)
  )])

  genealogy_new_med <- unique(data.table::rbindlist(
    list(
      top_parents_table,
      top_medium_table,
      low_medium_table_2,
      low_medium_table_1
    ),
    use.names = TRUE,
    fill = TRUE
  ))

  table_next <- data.table::copy(genealogy[
    !labels %in% genealogy_new_med$labels &
      parents %in% top_medium_table$ids &
      parents != ""
  ])
  if (nrow(table_next) > 0L) {
    table_next[,
      new_labels := trimws(paste(
        gsub("^([^-]*)-(.*)", "", parents),
        labels
      ))
    ]
    table_next[,
      ids := paste(
        gsub("^([^-]*)-(.*)", "\\2", parents),
        new_labels,
        sep = "-"
      )
    ]
    table_next[grepl("^Other", parents), new_labels := labels]
  } else {
    table_next[, new_labels := character(0)]
  }

  genealogy_new_med_2 <- unique(data.table::rbindlist(
    list(genealogy_new_med, table_next),
    use.names = TRUE,
    fill = TRUE
  ))

  tc_labels <- genealogy[!labels %in% genealogy_new_med$ids, list(labels)]
  tc_join1 <- dt_left_join(tc_labels, genealogy_new_med_2)
  tc_join1 <- tc_join1[grepl("-", parents)]
  tc_join2 <- dt_inner_join(
    tc_join1,
    genealogy,
    by = c("parents" = "ids")
  )
  table_children <- tc_join2[, list(
    chemical_pathway = chemical_pathway.x,
    parents = parents,
    ids = ids,
    labels = labels.x,
    new_labels = new_labels
  )]

  table_other <- unique(genealogy_new_med_2[
    grepl(" Other-", ids),
    list(chemical_pathway, ids = parents)
  ])
  if (nrow(table_other) > 0L) {
    table_other[, `:=`(
      parents = gsub("-.*", "\\1", ids),
      labels = gsub(".*-", "\\1", ids)
    )]
    table_other[, new_labels := labels]
  } else {
    table_other[, `:=`(
      parents = character(0),
      labels = character(0),
      new_labels = character(0)
    )]
  }

  genealogy_new_med_3 <- unique(data.table::rbindlist(
    list(genealogy_new_med_2, table_children, table_other),
    use.names = TRUE,
    fill = TRUE
  ))

  df2_23 <- unique(dataframe2[, list(
    best_candidate_3 = as.character(best_candidate_3),
    best_candidate_2 = as.character(best_candidate_2)
  )])

  mc1_base <- unique(genealogy_new_med_3[
    grepl(" Other$", parents),
    list(chemical_pathway, parents = ids, best_candidate_2 = labels)
  ])
  mc1_join <- dt_left_join(mc1_base, df2_23)
  missing_children_1 <- mc1_join[, list(
    chemical_pathway = as.character(chemical_pathway),
    parents = parents,
    ids = paste(best_candidate_2, best_candidate_3, sep = "-"),
    labels = best_candidate_3,
    new_labels = best_candidate_3
  )]
  for (col in names(missing_children_1)) {
    data.table::set(
      missing_children_1,
      j = col,
      value = gsub("-NA$", "", missing_children_1[[col]])
    )
  }

  mc2_base <- unique(genealogy_new_med_3[
    grepl("^Other", parents),
    list(chemical_pathway, parents = ids, best_candidate_2 = labels)
  ])
  mc2_join <- dt_left_join(mc2_base, df2_23)
  missing_children_2 <- unique(mc2_join[, list(
    chemical_pathway = as.character(chemical_pathway),
    parents = parents,
    ids = paste(gsub(".*-", "\\1", parents), best_candidate_3, sep = "-"),
    labels = best_candidate_3,
    new_labels = best_candidate_3
  )])
  for (col in names(missing_children_2)) {
    data.table::set(
      missing_children_2,
      j = col,
      value = gsub("-NA$", "", missing_children_2[[col]])
    )
  }
  missing_children_2 <- missing_children_2[!is.na(labels)]

  missing_children <- data.table::rbindlist(
    list(missing_children_1, missing_children_2),
    use.names = TRUE,
    fill = TRUE
  )

  g3_char <- data.table::copy(genealogy_new_med_3)
  g3_char[, chemical_pathway := as.character(chemical_pathway)]
  genealogy_new_med_4 <- unique(data.table::rbindlist(
    list(g3_char, missing_children),
    use.names = TRUE,
    fill = TRUE
  )[, list(chemical_pathway, parents, ids, labels, new_labels)])
  genealogy_new_med_4[, link := gsub("-.*", "", ids)]

  df2_sp <- data.table::copy(dataframe2[!is.na(species)])
  df2_sp[, `:=`(
    best_candidate_3 = as.character(best_candidate_3),
    chemical_pathway = as.character(chemical_pathway)
  )]

  table_new <- dt_full_join(
    df2_sp,
    genealogy_new_med_4,
    by = c(
      "best_candidate_2" = "link",
      "best_candidate_3" = "labels",
      "chemical_pathway" = "chemical_pathway"
    )
  )

  if (type == "analysis") {
    if (nrow(table_new) > 0L) {
      table_new <- table_new[is.na(id) | (!is.na(sample) & id == sample)]
    } else {
      table_new[, `:=`(integral = NA_real_, new_labels = NA_character_)]
    }
    if (detector == "cad") {
      table_new[, intensity := integral]
    }
    table_new <- unique(table_new[, list(
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
    )])
    table_new[, values := .N, by = c("parents", "ids", "sample", "intensity")]
    table_1_new <- unique(table_new[
      !is.na(ids),
      list(parents, ids, labels, values, sample, intensity, species)
    ])
  } else {
    table_new[, labels := best_candidate_3]
    table_new[, values := .N, by = c("parents", "ids", "sample")]
    table_1_new <- unique(table_new[
      !is.na(ids),
      list(parents, ids, labels, values, sample, organism, species)
    ])
  }

  table_1_new[, labels := as.character(labels)]
  additional_row <- data.table::copy(table_1_new[labels == "Other other"])
  if (nrow(additional_row) > 0L) {
    additional_row[, `:=`(
      parents = "Other-Other other",
      ids = "Other other - Other other other",
      labels = "Other other other"
    )]
  }
  table_1_1_new <- data.table::rbindlist(
    list(table_1_new, additional_row),
    use.names = TRUE,
    fill = TRUE
  )

  final_table_3_1 <- unique(table_1_1_new[
    ids %in% low_medium_table_1_a$ids,
    list(parents, ids, labels)
  ])

  val_col <- if (type == "analysis") "intensity" else "values"
  ft4_sub <- table_1_1_new[
    grepl("-", parents) & parents %in% missing_children_1$parents
  ]
  if (nrow(ft4_sub) > 0L) {
    final_table_4_1 <- ft4_sub[,
      list(values = sum(get(val_col))),
      keyby = c("parents", "ids", "labels", "sample", "species")
    ][!is.na(labels)]
    data.table::setkey(final_table_4_1, NULL)
    final_table_4_1[, join := gsub("-.*", "", parents)]
  } else {
    final_table_4_1 <- data.table::data.table(
      parents = character(0),
      ids = character(0),
      labels = character(0),
      sample = character(0),
      species = character(0),
      values = numeric(0),
      join = character(0)
    )
  }

  final_table_4 <- final_table_4_1[, list(
    parents,
    ids,
    labels,
    sample,
    species,
    values
  )]

  ft3_2_join <- dt_inner_join(
    final_table_3_1,
    final_table_4_1,
    by = c("labels" = "join")
  )
  if (nrow(ft3_2_join) > 0L) {
    ft3_2_join[,
      labels_clean := gsub(
        "([a-zA-Z]|\\))(-)([a-zA-Z]|[0-9])(.*)",
        "\\1",
        ids.y
      )
    ]
    final_table_3_2 <- ft3_2_join[,
      list(values = sum(values)),
      keyby = list(
        parents = ids.x,
        ids = parents.y,
        labels = labels_clean,
        sample = sample,
        species = species
      )
    ]
    data.table::setkey(final_table_3_2, NULL)
  } else {
    final_table_3_2 <- data.table::data.table(
      parents = character(0),
      ids = character(0),
      labels = character(0),
      sample = character(0),
      species = character(0),
      values = numeric(0)
    )
  }

  ft3_3_sub <- table_1_1_new[
    !is.na(species) &
      grepl("-", parents) &
      !parents %in% missing_children_1$parents
  ]
  if (nrow(ft3_3_sub) > 0L) {
    final_table_3_3 <- ft3_3_sub[,
      list(values = sum(get(val_col))),
      keyby = c("parents", "ids", "labels", "sample", "species")
    ]
    data.table::setkey(final_table_3_3, NULL)
  } else {
    final_table_3_3 <- data.table::data.table(
      parents = character(0),
      ids = character(0),
      labels = character(0),
      sample = character(0),
      species = character(0),
      values = numeric(0)
    )
  }

  final_table_3 <- unique(data.table::rbindlist(
    list(final_table_3_3, final_table_3_2),
    use.names = TRUE,
    fill = TRUE
  ))
  if (nrow(final_table_3) > 0L) {
    data.table::setorder(final_table_3, parents, ids, na.last = TRUE)
  }

  t11_base <- table_1_1_new[, list(parents, ids, labels)]
  ft2_join <- dt_left_join(
    t11_base,
    final_table_3[, list(values, sample, parents, species)],
    by = c("ids" = "parents")
  )
  final_table_2 <- ft2_join[,
    list(values = sum(values)),
    keyby = c("parents", "ids", "labels", "sample", "species")
  ][!is.na(values)]
  data.table::setkey(final_table_2, NULL)

  ft1_join <- dt_left_join(
    t11_base,
    final_table_2[, list(values, sample, parents, species)],
    by = c("ids" = "parents")
  )
  final_table_1 <- ft1_join[
    !is.na(values),
    list(values = sum(values)),
    keyby = c("parents", "ids", "labels", "sample", "species")
  ]
  data.table::setkey(final_table_1, NULL)

  final_table <- data.table::rbindlist(
    list(final_table_1, final_table_2, final_table_3, final_table_4),
    use.names = TRUE,
    fill = TRUE
  )[!is.na(sample)]

  if (type == "analysis" && nrow(final_table) > 0L) {
    if (rescale) {
      final_table[,
        values := values / (sum(values) / 3),
        by = c("sample", "species")
      ]
    } else {
      final_table[, values := values / 3, by = "sample"]
    }
  }

  final_table
}

#' Build literature hierarchies across a list of tables
#'
#' @param xs List of data frames / `data.table` objects.
#' @param comparison Optional comparison taxon string(s).
#'
#' @return A list of hierarchy `data.table` objects.
#' @keywords internal
#' @noRd
hierarchies_progress <- function(xs, comparison = NULL) {
  has_multiword_comp <- !is.null(comparison) &&
    !all(is.na(comparison)) &&
    all(grepl("\\W+", comparison))

  lapply(xs, function(x) {
    dt <- data.table::copy(data.table::as.data.table(x))
    if (nrow(dt) != 0L) {
      dt[, `:=`(
        best_candidate_1 = chemical_pathway,
        best_candidate_2 = chemical_superclass,
        best_candidate_3 = chemical_class,
        organism = if (has_multiword_comp) {
          taxaLabels
        } else {
          gsub(" .*", "", taxaLabels)
        }
      )]
      dt[, `:=`(sample = organism, species = organism)]
      drop_cols <- intersect(
        c("taxa", "taxaLabels", "references", "referencesLabels"),
        names(dt)
      )
      if (length(drop_cols) > 0L) {
        dt[, (drop_cols) := NULL]
      }
      prepare_hierarchy(
        dataframe = unique(dt),
        type = "literature"
      )
    } else {
      data.table::data.table(
        parents = NA_character_,
        ids = NA_character_,
        labels = NA_character_,
        sample = NA_character_,
        species = NA_character_,
        values = 0
      )
    }
  })
}

#' Build grouped literature hierarchies across a list of tables
#'
#' @param xs List of data frames / `data.table` objects.
#'
#' @return A list of grouped hierarchy `data.table` objects.
#' @keywords internal
#' @noRd
hierarchies_grouped_progress <- function(xs) {
  lapply(xs, function(x) {
    dt <- data.table::copy(data.table::as.data.table(x))
    if (nrow(dt) != 0L) {
      non_var_taxa <- unique(dt[!grepl("var. ", taxaLabels), list(taxa)])
      if (nrow(non_var_taxa) != 1L) {
        dt[, `:=`(
          best_candidate_1 = chemical_pathway,
          best_candidate_2 = chemical_superclass,
          best_candidate_3 = chemical_class,
          organism = taxaLabels
        )]
        dt[, `:=`(sample = organism, species = organism)]
        drop_cols <- intersect(
          c("taxa", "taxaLabels", "references", "referencesLabels"),
          names(dt)
        )
        if (length(drop_cols) > 0L) {
          dt[, (drop_cols) := NULL]
        }
        prepare_hierarchy(
          dataframe = unique(dt),
          type = "literature"
        )
      } else {
        NULL
      }
    } else {
      data.table::data.table(
        parents = NA_character_,
        ids = NA_character_,
        labels = NA_character_,
        sample = NA_character_,
        species = NA_character_,
        values = 0
      )
    }
  })
}

#' Prepare hierarchy table for stacked bar plots
#'
#' @param dataframe Hierarchy `data.table` from `prepare_hierarchy()`.
#' @param organism Column name representing the organism (default `"species"`).
#'
#' @return A `data.table` with factor-ordered `ids`, `color`, and `sample`
#'   columns ready for plotting.
#' @keywords internal
#' @noRd
prepare_plot <- function(dataframe, organism = "species") {
  dt <- data.table::copy(data.table::as.data.table(dataframe))
  presamples <- dt[
    parents != "" &
      !grepl("-", parents) &
      !is.na(get(organism))
  ]
  if (nrow(presamples) == 0L) {
    return(presamples)
  }

  presamples[, species := get(organism)]
  presamples[, valuez := sum(values), by = "parents"]
  data.table::setorder(presamples, -valuez, -values, na.last = TRUE)
  presamples[, group := as.integer(factor(parents, levels = unique(parents)))]
  presamples[,
    subgroup := dense_rank_base(as.integer(factor(ids, levels = unique(ids)))),
    by = "group"
  ]
  presamples[, tot := sum(values), by = "sample"]

  tempval <- unique(presamples[parents == "Other", group])
  if (length(tempval) != 0L) {
    tv <- tempval[1L]
    samples_0 <- data.table::copy(presamples[group == tv])
    samples_0[, group := 666]
    samples_1 <- data.table::copy(presamples[group > tv])
    samples_1[, group := group - 1]
    samples_2 <- presamples[group < tv]
    samples <- data.table::rbindlist(list(samples_0, samples_1, samples_2))
  } else {
    samples <- presamples
  }

  is_other <- grepl("Other", samples$parents, fixed = TRUE)
  colors <- character(nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sg <- samples$subgroup[i]
    if (is_other[i]) {
      pal <- microshades_grey[[1L]]
      colors[i] <- if (!is.na(sg) && sg >= 1L && sg <= length(pal)) {
        pal[[sg]]
      } else {
        NA_character_
      }
    } else {
      grp <- samples$group[i]
      if (!is.na(grp) && grp >= 1L && grp <= length(microshades)) {
        pal <- microshades[[grp]]
        colors[i] <- if (!is.na(sg) && sg >= 1L && sg <= length(pal)) {
          pal[[sg]]
        } else {
          NA_character_
        }
      } else {
        colors[i] <- NA_character_
      }
    }
  }
  samples[, `:=`(
    color = colors,
    relative = values / tot
  )]
  data.table::setorder(samples, group, subgroup, na.last = TRUE)
  samples[, idz := paste(group, subgroup)]

  na_col <- is.na(samples$color)
  if (any(na_col)) {
    quickfix <- data.table::copy(samples[na_col])
    quickfix[, subsub := .GRP, by = c("group", "subgroup")]
    quickfix[, color := paste0("#00", sprintf("%04d", subsub))]
    quickfix[, subsub := NULL]
    samples <- data.table::rbindlist(list(samples[!na_col], quickfix))
  }

  samples[, ids := fct_reorder2_base(ids, group, idz, desc = FALSE)]
  samples[, color := fct_reorder2_base(color, group, idz, desc = FALSE)]
  samples[, sample := fct_reorder2_base(sample, valuez, sample, desc = TRUE)]

  samples
}

#' Prepare feature comparison table for pseudochromatogram bar plots
#'
#' @param dataframe Data frame or `data.table` of features.
#'
#' @return A `data.table` with chemical class colors and biological source
#'   levels ready for `plot_histograms_confident()` and `plot_histograms_taxo()`.
#' @keywords internal
#' @noRd
prepare_plot_2 <- function(dataframe) {
  dt <- data.table::copy(data.table::as.data.table(dataframe))
  dt <- dt[!is.na(inchikey_2D)]
  dt <- unique(
    dt,
    by = c(
      "peak_id",
      "inchikey_2D",
      "best_candidate_1",
      "best_candidate_2",
      "best_candidate_3"
    )
  )
  if (nrow(dt) == 0L) {
    return(dt)
  }

  dt[, peak_area := peak_area / .N, by = "peak_id"]
  dt[
    is.na(best_candidate_1) |
      best_candidate_1 %in% c("notConfident", "notClassified"),
    best_candidate_1 := "Other"
  ]
  dt[is.na(best_candidate_2), best_candidate_2 := "notAnnotated"]
  dt[
    best_candidate_2 == "notConfident notConfident",
    best_candidate_2 := "notConfident"
  ]
  dt[is.na(best_candidate_3), best_candidate_3 := "notAnnotated notAnnotated"]
  dt[
    best_candidate_3 == "notConfident notConfident notConfident",
    best_candidate_3 := "notConfident notConfident"
  ]
  dt[
    best_candidate_3 == "notClassified",
    best_candidate_3 := "notClassified notClassified"
  ]

  dt[, sum := sum(unique(peak_area)), by = "best_candidate_1"]
  data.table::setorder(dt, -sum, na.last = TRUE)
  dt[,
    group := as.integer(factor(
      best_candidate_1,
      levels = unique(best_candidate_1)
    ))
  ]
  dt[,
    subgroup := dense_rank_base(as.integer(factor(
      best_candidate_2,
      levels = unique(best_candidate_2)
    ))),
    by = "group"
  ]

  is_other <- grepl("Other", dt$best_candidate_1, fixed = TRUE)
  colors <- character(nrow(dt))
  for (i in seq_len(nrow(dt))) {
    sg <- dt$subgroup[i]
    if (is_other[i]) {
      colors[i] <- microshades_grey[[1L]][[sg]]
    } else {
      colors[i] <- microshades[[dt$group[i]]][[sg]]
    }
  }
  dt[, `:=`(
    color = colors,
    name = paste(best_candidate_1, best_candidate_2, sep = " - ")
  )]

  dataframe_prep <- dt[, list(
    sample = sample,
    sample_organism = species,
    peak_id = peak_id,
    peak_rt = peak_rt_apex,
    peak_area = peak_area,
    feature_id = feature_id,
    feature_area = feature_area,
    feature_rt = rt,
    feature_mz = mz,
    comparison_score = comparison_score,
    smiles_2D = smiles_2D,
    inchikey_2D = inchikey_2D,
    score_biological = score_biological,
    score_chemical = score_chemical,
    score_final = score_final,
    nearest_organism = best_candidate_organism,
    best_candidate_1 = best_candidate_1,
    best_candidate_2 = best_candidate_2,
    best_candidate_3 = best_candidate_3,
    group = group,
    subgroup = subgroup,
    color = color,
    name = name,
    sum = sum
  )]

  sb <- suppressWarnings(as.numeric(dataframe_prep$score_biological))
  name_2_val <- ifelse(
    !is.na(sb) & sb >= 0.9,
    "Species",
    ifelse(
      !is.na(sb) & sb >= 0.8,
      "Genus",
      ifelse(
        !is.na(sb) & sb >= 0.6,
        "Family",
        ifelse(!is.na(sb) & sb >= 0.2, "Kingdom", "Other")
      )
    )
  )
  color_2_val <- ifelse(
    name_2_val == "Species",
    microshades[[9L]][[1L]],
    ifelse(
      name_2_val == "Genus",
      microshades[[9L]][[2L]],
      ifelse(
        name_2_val == "Family",
        microshades[[9L]][[3L]],
        ifelse(
          name_2_val == "Kingdom",
          microshades[[9L]][[4L]],
          microshades_grey[[1L]][[2L]]
        )
      )
    )
  )

  dataframe_prep[, `:=`(
    name_2 = name_2_val,
    color_2 = color_2_val
  )]

  if (nrow(dataframe_prep) > 0L) {
    dataframe_prep[, `:=`(
      color = fct_reorder2_base(color, sum, group, desc = TRUE),
      name = fct_reorder2_base(name, sum, group, desc = TRUE),
      color_2 = fct_reorder2_base(
        color_2,
        score_biological,
        score_biological,
        desc = TRUE
      ),
      name_2 = fct_reorder2_base(
        name_2,
        score_biological,
        score_biological,
        desc = TRUE
      )
    )]
  }

  dataframe_prep
}

#' Plot pseudochromatogram colored by chemical pathway and superclass
#'
#' @param dataframe Prepared `data.table` from `prepare_plot_2()`.
#' @param chromatogram Chromatogram `data.table` with `rtime` and `intensity`.
#' @param level Either `"max"` (peak level) or `"min"` (feature level).
#' @param time_min Minimum retention time (minutes).
#' @param time_max Maximum retention time (minutes).
#'
#' @return A `ggplot` object.
#' @keywords internal
#' @noRd
plot_histograms_confident <- function(
  dataframe,
  chromatogram,
  level = "max",
  time_min,
  time_max
) {
  dt <- data.table::copy(data.table::as.data.table(dataframe))
  dt[, n := .N, by = "peak_area"]
  dt[, m := .N, by = "feature_area"]

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = chromatogram,
      mapping = ggplot2::aes(x = rtime, y = intensity / max(intensity)),
      col = "black",
      linewidth = 0.1
    ) +
    ggplot2::geom_bar(
      data = dt,
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
    ) +
    ggplot2::scale_fill_manual(
      values = as.character(levels(dt$color)),
      guide = ggplot2::guide_legend(ncol = 1)
    ) +
    ggplot2::labs(fill = "Chemical Pathway - Superclass") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text.y = ggplot2::element_text(face = "italic")
    ) +
    ggplot2::ylab("Intensity") +
    ggplot2::xlab("Retention time [min]") +
    ggplot2::xlim(max(time_min), min(time_max))
}

#' Plot pseudochromatogram colored by closest taxonomic occurrence
#'
#' @param dataframe Prepared `data.table` from `prepare_plot_2()`.
#' @param chromatogram Chromatogram `data.table` with `rtime` and `intensity`.
#' @param level Either `"max"` (peak level) or `"min"` (feature level).
#' @param mode Ionization mode (`"pos"` or `"neg"`).
#' @param time_min Minimum retention time (minutes).
#' @param time_max Maximum retention time (minutes).
#'
#' @return A `ggplot` object.
#' @keywords internal
#' @noRd
plot_histograms_taxo <- function(
  dataframe,
  chromatogram,
  level = "max",
  mode = "pos",
  time_min,
  time_max
) {
  dt <- data.table::copy(data.table::as.data.table(dataframe))
  chrom <- data.table::copy(data.table::as.data.table(chromatogram))
  dt[, n := .N, by = "peak_area"]
  dt[, m := .N, by = "feature_area"]

  if (mode == "neg") {
    dt[, `:=`(
      peak_area = -1 * peak_area,
      feature_area = -1 * feature_area
    )]
    chrom[, intensity := -1 * intensity]
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = chrom,
      mapping = ggplot2::aes(x = rtime, y = intensity / max(abs(intensity))),
      col = "black",
      linewidth = 0.1
    ) +
    ggplot2::geom_bar(
      data = dt,
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
    ) +
    ggplot2::scale_fill_manual(
      values = as.character(levels(dt$color_2))
    ) +
    ggplot2::labs(fill = "Already reported in") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text.y = ggplot2::element_text(face = "italic")
    ) +
    ggplot2::ylab("Intensity") +
    ggplot2::xlab("Retention time [min]") +
    ggplot2::xlim(max(time_min), min(time_max))
}

#' Plot horizontal stacked bar chart of chemical classes across samples/taxa
#'
#' @param dataframe Prepared `data.table` from `prepare_plot()`.
#' @param label X-axis label.
#' @param y Column name for bar heights (default `"values"`).
#' @param xlab Logical; whether to display discrete X-axis labels.
#'
#' @return A `ggplot` object.
#' @keywords internal
#' @noRd
plot_histograms_litt <- function(
  dataframe,
  label,
  y = "values",
  xlab = TRUE
) {
  ggplot2::ggplot(
    dataframe,
    ggplot2::aes(x = sample, y = get(y), fill = ids)
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::scale_fill_manual(
      values = as.character(levels(dataframe$color)),
      guide = ggplot2::guide_legend(ncol = 1)
    ) +
    {
      if (xlab) {
        ggplot2::scale_x_discrete(labels = levels(dataframe$sample))
      }
    } +
    ggplot2::labs(fill = "Chemical class") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(face = "italic")
    ) +
    ggplot2::xlab(label) +
    ggplot2::ylab("Count") +
    ggplot2::coord_flip()
}

#' Generate pseudochromatogram plot suite across taxonomic and confidence levels
#'
#' @param list Output list from `prepare_comparison()`.
#' @param chromatogram Baseline-corrected chromatogram `data.table`.
#' @param mode Ionization mode (`"pos"` or `"neg"`).
#' @param time_min Minimum retention time (minutes).
#' @param time_max Maximum retention time (minutes).
#'
#' @return Named list of 8 `ggplot` pseudochromatograms.
#' @keywords internal
#' @noRd
plot_results_1 <- function(
  list,
  chromatogram,
  mode = "pos",
  time_min,
  time_max
) {
  message("preparing histograms")
  maj_other <- make_other(list$peaks_maj_precor_taxo_cor)
  min_other <- make_other(list$peaks_min_precor_taxo_cor)

  df_histogram_maj <- prepare_plot_2(maj_other)
  df_histogram_min <- prepare_plot_2(min_other)
  df_histogram_maj_conf <- prepare_plot_2(no_other(maj_other))
  df_histogram_min_conf <- prepare_plot_2(no_other(min_other))

  message("plotting histograms...")
  message("... taxo")
  histograms_taxo_maj <- plot_histograms_taxo(
    df_histogram_maj,
    chromatogram = chromatogram,
    mode = mode,
    time_min = time_min,
    time_max = time_max
  )
  histograms_taxo_min <- plot_histograms_taxo(
    df_histogram_min,
    level = "min",
    chromatogram = chromatogram,
    mode = mode,
    time_min = time_min,
    time_max = time_max
  )

  message("... confident features")
  histograms_conf_maj <- plot_histograms_confident(
    df_histogram_maj_conf,
    chromatogram = chromatogram,
    time_min = time_min,
    time_max = time_max
  )
  histograms_conf_min <- plot_histograms_confident(
    df_histogram_min_conf,
    level = "min",
    chromatogram = chromatogram,
    time_min = time_min,
    time_max = time_max
  )

  message("... unique structures")
  histograms_unique_maj <- plot_histograms_confident(
    unique(df_histogram_maj, by = c("peak_id", "inchikey_2D")),
    chromatogram = chromatogram,
    time_min = time_min,
    time_max = time_max
  )
  histograms_unique_min <- plot_histograms_confident(
    unique(df_histogram_min, by = c("peak_id", "inchikey_2D")),
    level = "min",
    chromatogram = chromatogram,
    time_min = time_min,
    time_max = time_max
  )

  message("... confident structures")
  histograms_unique_conf_maj <- plot_histograms_confident(
    unique(df_histogram_maj_conf, by = c("peak_id", "inchikey_2D")),
    chromatogram = chromatogram,
    time_min = time_min,
    time_max = time_max
  )
  histograms_unique_conf_min <- plot_histograms_confident(
    unique(df_histogram_min_conf, by = c("peak_id", "inchikey_2D")),
    level = "min",
    chromatogram = chromatogram,
    time_min = time_min,
    time_max = time_max
  )

  list(
    histograms_taxo_maj = histograms_taxo_maj,
    histograms_taxo_min = histograms_taxo_min,
    histograms_conf_maj = histograms_conf_maj,
    histograms_conf_min = histograms_conf_min,
    histograms_unique_maj = histograms_unique_maj,
    histograms_unique_min = histograms_unique_min,
    histograms_unique_conf_maj = histograms_unique_conf_maj,
    histograms_unique_conf_min = histograms_unique_conf_min
  )
}

#' Generate signal-based sunburst plots from peak comparison tables
#'
#' @param list Output list from `prepare_comparison()`.
#'
#' @return A named list of two `plotly` sunburst widgets:
#'   `sunburst_signal_based_duo` and `sunburst_conf_signal_based_duo`.
#' @keywords internal
#' @noRd
plot_results_2 <- function(list) {
  dedup_peaks <- function(df, colname = "peak_rt_apex") {
    dt <- data.table::copy(data.table::as.data.table(df))
    if (nrow(dt) == 0L) {
      return(dt)
    }
    data.table::setorder(dt, -comparison_score, na.last = TRUE)
    unique(
      dt,
      by = c(
        colname,
        "inchikey_2D",
        "best_candidate_1",
        "best_candidate_2",
        "best_candidate_3"
      )
    )
  }

  strip_posneg <- function(df) {
    dt <- data.table::copy(data.table::as.data.table(df))
    if (nrow(dt) == 0L) {
      return(dt)
    }
    dt[, `:=`(
      newrt = round(peak_rt_apex, 1),
      sample = gsub("_(pos|neg)", "", sample, ignore.case = TRUE),
      id = gsub("_(pos|neg)", "", id, ignore.case = TRUE)
    )]
    dt
  }

  quick_sb <- function(df) {
    plotly::plot_ly(
      data = df,
      ids = ~ids,
      labels = ~labels,
      parents = ~parents,
      values = ~values,
      maxdepth = 3,
      type = "sunburst",
      branchvalues = "total",
      textinfo = "label+percent value+percent parent+percent root"
    )
  }

  message("preparing hierarchies...")
  message("... on everything")
  table_taxo_maj_cor_signal_2 <- prepare_hierarchy(
    dedup_peaks(
      strip_posneg(list$peaks_maj_precor_taxo_cor),
      colname = "newrt"
    ),
    detector = "cad"
  )

  message("... on confident")
  table_taxo_maj_cor_conf_signal_2 <- prepare_hierarchy(
    dedup_peaks(
      strip_posneg(no_other(list$peaks_maj_precor_taxo_cor)),
      colname = "newrt"
    ),
    detector = "cad"
  )

  index_signal_duo <- table_taxo_maj_cor_signal_2[parents == ""]
  if (nrow(index_signal_duo) > 0L) {
    data.table::setorder(index_signal_duo, -values, na.last = TRUE)
    index_signal_duo <- unique(index_signal_duo, by = c("ids", "labels"))
  }

  sunburst_microshades_grey_signal_duo <- microshades_colors
  index_s_duo <- which(index_signal_duo$ids == "Other")
  if (length(index_s_duo) > 0L) {
    sunburst_microshades_grey_signal_duo[index_s_duo] <- microshades_grey[[
      1L
    ]][[5L]]
  }

  sunburst_signal_based_duo <- plotly::layout(
    quick_sb(table_taxo_maj_cor_signal_2),
    colorway = sunburst_microshades_grey_signal_duo
  )

  sunburst_conf_signal_based_duo <- plotly::layout(
    quick_sb(table_taxo_maj_cor_conf_signal_2),
    colorway = microshades_colors
  )

  list(
    sunburst_signal_based_duo = sunburst_signal_based_duo,
    sunburst_conf_signal_based_duo = sunburst_conf_signal_based_duo
  )
}

#' Build interactive plotly treemaps with structure count titles
#'
#' @param xs Character vector of hierarchy names.
#' @param type Plotly trace type (default `"treemap"`).
#' @param hierarchies Named list of hierarchy `data.table` objects.
#' @param tables Named list of compound-taxon `data.table` objects used to
#'   count distinct structures for the plot title.
#'
#' @return A named list of titled `plotly` treemap/sunburst widgets.
#' @keywords internal
#' @noRd
treemaps_progress <- function(xs, type = "treemap", hierarchies, tables) {
  count_structures <- function(nm) {
    if (is.null(tables) || is.null(tables[[nm]])) {
      return(0L)
    }
    t_dt <- data.table::as.data.table(tables[[nm]])
    if (!"structure" %in% names(t_dt)) {
      return(nrow(t_dt))
    }
    nrow(unique(t_dt[, list(structure)]))
  }

  lapply(stats::setNames(xs, xs), function(x) {
    h_dt <- data.table::as.data.table(hierarchies[[x]])
    if (x == "special") {
      sp_vals <- unique(h_dt$species)
      d1 <- data.table::as.data.table(hierarchies[[sp_vals[1L]]])[
        labels != "notClassified sub"
      ]
      d2 <- data.table::as.data.table(hierarchies[[sp_vals[2L]]])[
        labels != "notClassified sub"
      ]
      plotly::plot_ly() |>
        plotly::add_trace(
          data = d1,
          ids = ~ids,
          labels = ~labels,
          parents = ~parents,
          values = ~values,
          maxdepth = 3,
          type = type,
          branchvalues = "total",
          textinfo = "label+percent value+percent parent+percent root",
          domain = list(row = 0, column = 0)
        ) |>
        plotly::add_trace(
          data = d2,
          ids = ~ids,
          labels = ~labels,
          parents = ~parents,
          values = ~values,
          maxdepth = 3,
          type = type,
          branchvalues = "total",
          textinfo = "label+percent value+percent parent+percent root",
          domain = list(row = 0, column = 1)
        ) |>
        plotly::layout(
          title = paste(
            "Comparative analysis",
            "\n",
            sp_vals[1L],
            "(",
            count_structures(sp_vals[1L]),
            ")",
            "                                 ",
            sp_vals[2L],
            "(",
            count_structures(sp_vals[2L]),
            ")"
          ),
          grid = list(rows = 1, columns = 2),
          colorway = microshades_colors,
          margin = list(t = 40)
        )
    } else {
      plotly::plot_ly(
        data = h_dt[labels != "notClassified sub"],
        ids = ~ids,
        labels = ~labels,
        parents = ~parents,
        values = ~values,
        maxdepth = 3,
        type = type,
        branchvalues = "total",
        textinfo = "label+percent value+percent parent+percent root"
      ) |>
        plotly::layout(
          colorway = microshades_colors,
          title = paste(x, "(", count_structures(x), ")"),
          margin = list(t = 40)
        )
    }
  })
}

#' Build interactive plotly treemaps from hierarchy tables
#'
#' @param xs Character vector of hierarchy names.
#' @param type Plotly trace type (default `"treemap"`).
#' @param hierarchies Named list of hierarchy `data.table` objects.
#'
#' @return A named list of `plotly` treemap widgets.
#' @keywords internal
#' @noRd
treemaps_progress_no_title <- function(xs, type = "treemap", hierarchies) {
  lapply(stats::setNames(xs, xs), function(x) {
    h_dt <- data.table::as.data.table(hierarchies[[x]])
    if (x == "special") {
      sp_vals <- unique(h_dt$species)
      d1 <- data.table::as.data.table(hierarchies[[sp_vals[1L]]])[
        labels != "notClassified sub"
      ]
      d2 <- data.table::as.data.table(hierarchies[[sp_vals[2L]]])[
        labels != "notClassified sub"
      ]
      plotly::plot_ly() |>
        plotly::add_trace(
          data = d1,
          ids = ~ids,
          labels = ~labels,
          parents = ~parents,
          values = ~values,
          maxdepth = 3,
          type = type,
          textinfo = "label+percent value+percent parent+percent root",
          domain = list(row = 0, column = 0)
        ) |>
        plotly::add_trace(
          data = d2,
          ids = ~ids,
          labels = ~labels,
          parents = ~parents,
          values = ~values,
          maxdepth = 3,
          type = type,
          textinfo = "label+percent value+percent parent+percent root",
          domain = list(row = 0, column = 1)
        ) |>
        plotly::layout(
          grid = list(rows = 1, columns = 2),
          colorway = microshades_colors,
          margin = list(t = 40)
        )
    } else {
      plotly::plot_ly(
        data = h_dt[labels != "notClassified sub"],
        ids = ~ids,
        labels = ~labels,
        parents = ~parents,
        values = ~values,
        maxdepth = 3,
        type = type,
        branchvalues = "total",
        textinfo = "label+percent value+percent parent+percent root"
      ) |>
        plotly::layout(
          colorway = microshades_colors,
          margin = list(t = 40)
        )
    }
  })
}
