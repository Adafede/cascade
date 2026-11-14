#' Replace matching values with NA
#'
#' @param x A vector or factor.
#' @param y Value to replace with `NA`.
#'
#' @return Vector `x` with occurrences of `y` replaced by `NA`.
#' @keywords internal
#' @noRd
y_as_na <- function(x, y) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  x[as.character(x) == y] <- NA
  x
}

#' URL-encode unique SMILES strings once and map back
#'
#' @param x Character vector of SMILES strings.
#' @param empty Character(1) placeholder returned for missing/empty inputs.
#' @param reserved Logical(1) passed to [utils::URLencode()].
#'
#' @return Character vector of URL-encoded SMILES strings.
#' @keywords internal
#' @noRd
encode_smiles <- function(x, empty = NA_character_, reserved = TRUE) {
  x_chr <- as.character(x)
  out <- rep(empty, length(x_chr))
  valid <- !is.na(x_chr) & nzchar(x_chr) & x_chr != "NA"
  if (any(valid)) {
    uniq_s <- unique(x_chr[valid])
    enc_s <- vapply(
      uniq_s,
      utils::URLencode,
      character(1L),
      reserved = reserved,
      USE.NAMES = FALSE
    )
    out[valid] <- enc_s[match(x_chr[valid], uniq_s)]
  }
  out
}

#' Order-preserving left join on data.tables (tidytable-compatible)
#'
#' @param x,y `data.table` or `data.frame` objects to join.
#' @param by Character vector of join column names (named or unnamed), or `NULL`
#'   for natural join.
#' @param suffix Character(2) suffixes appended to colliding non-key columns.
#'
#' @return A `data.table` preserving `x`'s row order and left-to-right column
#'   order.
#' @keywords internal
#' @noRd
dt_left_join <- function(x, y, by = NULL, suffix = c(".x", ".y")) {
  x <- data.table::copy(data.table::as.data.table(x))
  y <- data.table::copy(data.table::as.data.table(y))
  if (is.null(by)) {
    by_x <- by_y <- intersect(names(x), names(y))
  } else if (!is.null(names(by))) {
    by_x <- ifelse(names(by) == "", unname(by), names(by))
    by_y <- unname(by)
  } else {
    by_x <- by_y <- unname(by)
  }
  y_non <- setdiff(names(y), by_y)
  coll <- intersect(setdiff(names(x), by_x), y_non)
  if (length(coll) > 0L) {
    data.table::setnames(x, coll, paste0(coll, suffix[1L]))
    data.table::setnames(y, coll, paste0(coll, suffix[2L]))
  }
  data.table::setnames(y, by_y, by_x)
  out <- y[x, on = by_x, allow.cartesian = TRUE]
  data.table::setcolorder(out, c(names(x), setdiff(names(y), by_x)))
  out
}

#' Order-preserving inner join on data.tables (tidytable-compatible)
#'
#' @param x,y `data.table` or `data.frame` objects to join.
#' @param by Character vector of join column names (named or unnamed), or `NULL`
#'   for natural join.
#' @param suffix Character(2) suffixes appended to colliding non-key columns.
#'
#' @return A `data.table` preserving `y` lookup order and left-to-right column
#'   order.
#' @keywords internal
#' @noRd
dt_inner_join <- function(x, y, by = NULL, suffix = c(".x", ".y")) {
  x <- data.table::copy(data.table::as.data.table(x))
  y <- data.table::copy(data.table::as.data.table(y))
  if (is.null(by)) {
    by_x <- by_y <- intersect(names(x), names(y))
  } else if (!is.null(names(by))) {
    by_x <- ifelse(names(by) == "", unname(by), names(by))
    by_y <- unname(by)
  } else {
    by_x <- by_y <- unname(by)
  }
  y_non <- setdiff(names(y), by_y)
  coll <- intersect(setdiff(names(x), by_x), y_non)
  if (length(coll) > 0L) {
    data.table::setnames(x, coll, paste0(coll, suffix[1L]))
    data.table::setnames(y, coll, paste0(coll, suffix[2L]))
  }
  on_vec <- stats::setNames(by_y, by_x)
  x[y, on = on_vec, allow.cartesian = TRUE, nomatch = 0L]
}

#' Order-preserving full join on data.tables (tidytable-compatible)
#'
#' @param x,y `data.table` or `data.frame` objects to join.
#' @param by Character vector of join column names (named or unnamed), or `NULL`
#'   for natural join.
#' @param suffix Character(2) suffixes appended to colliding non-key columns.
#'
#' @return A `data.table` containing all matched and unmatched rows from `x` and
#'   `y` without key sorting.
#' @keywords internal
#' @noRd
dt_full_join <- function(x, y, by = NULL, suffix = c(".x", ".y")) {
  x <- data.table::copy(data.table::as.data.table(x))
  y <- data.table::copy(data.table::as.data.table(y))
  if (is.null(by)) {
    by_x <- by_y <- intersect(names(x), names(y))
  } else if (!is.null(names(by))) {
    by_x <- ifelse(names(by) == "", unname(by), names(by))
    by_y <- unname(by)
  } else {
    by_x <- by_y <- unname(by)
  }
  x_names <- names(x)
  y_non <- setdiff(names(y), by_y)
  coll <- intersect(setdiff(x_names, by_x), y_non)
  if (length(coll) > 0L) {
    x_names[match(coll, x_names)] <- paste0(coll, suffix[1L])
    y_non[match(coll, y_non)] <- paste0(coll, suffix[2L])
  }
  out <- merge(
    x = x,
    y = y,
    by.x = by_x,
    by.y = by_y,
    suffixes = suffix,
    all.x = TRUE,
    all.y = TRUE,
    allow.cartesian = TRUE,
    sort = FALSE
  )
  data.table::setkey(out, NULL)
  data.table::setcolorder(out, c(x_names, y_non))
  out
}

#' Evaluate quadratic acetonitrile response polynomial
#'
#' @param acn_eluent Numeric vector of acetonitrile percentage in eluent.
#' @param q1 Quadratic coefficient.
#' @param q2 Linear coefficient.
#' @param q3 Intercept coefficient.
#'
#' @return Numeric vector of polynomial values.
#' @keywords internal
#' @noRd
p_acn_i <- function(acn_eluent, q1, q2, q3) {
  q1 * acn_eluent^2 + q2 * acn_eluent + q3
}

#' Predict acetonitrile-corrected CAD response from peak area
#'
#' Applies the inverse quadratic calibration model to correct Charged Aerosol
#' Detector (CAD) peak areas for the mobile-phase acetonitrile percentage at
#' elution.
#'
#' @param acn Numeric vector of acetonitrile percentages (default `100`).
#' @param peak_area Numeric vector of raw peak areas.
#' @param p1q1,p1q2,p1q3 Quadratic coefficients for parameter 1.
#' @param p2q1,p2q2,p2q3 Quadratic coefficients for parameter 2.
#' @param p3q1,p3q2,p3q3 Quadratic coefficients for parameter 3.
#'
#' @return Numeric vector of corrected responses (concentrations).
#' @keywords internal
#' @noRd
predict_response <- function(
  acn = 100,
  peak_area,
  p1q1 = 0.00001,
  p1q2 = -0.0006,
  p1q3 = -0.0778,
  p2q1 = 0.00002,
  p2q2 = -0.00022,
  p2q3 = 0.05499,
  p3q1 = -0.00017,
  p3q2 = 0.0209,
  p3q3 = 1.4041
) {
  p_acn_1 <- p_acn_i(acn_eluent = acn, q1 = p1q1, q2 = p1q2, q3 = p1q3)
  p_acn_2 <- p_acn_i(acn_eluent = acn, q1 = p2q1, q2 = p2q2, q3 = p2q3)
  p_acn_3 <- p_acn_i(acn_eluent = acn, q1 = p3q1, q2 = p3q2, q3 = p3q3)
  denom <- 2 * p_acn_1 + p_acn_2

  ifelse(
    !is.na(peak_area) & peak_area > 0,
    peak_area^(1 / denom) * exp(-p_acn_3 / denom),
    0
  )
}

#' Format SMILES as CDK Depict web image for gt tables
#'
#' @param x Character vector of URL-encoded SMILES strings.
#'
#' @return HTML image tag generated by `gt::web_image()`.
#' @keywords internal
#' @noRd
molinfo <- function(x) {
  gt::web_image(
    url = paste0(
      "https://www.simolecule.com/cdkdepict/depict/bot/svg?smi=",
      x
    ),
    height = 75
  )
}

#' Keep best structural candidates from TIMA annotation table
#'
#' @param df Data frame or `data.table` of annotations.
#'
#' @return A `data.table` containing the top candidate per feature with
#'   standardized chemical class columns.
#' @keywords internal
#' @noRd
keep_best_candidates <- function(df) {
  dt <- data.table::copy(data.table::as.data.table(df))
  for (col in names(dt)) {
    val <- dt[[col]]
    if (!is.character(val)) {
      val <- as.character(val)
    }
    pos <- regexpr("|", val, fixed = TRUE)
    hit <- which(pos > 0L)
    if (length(hit) > 0L) {
      val[hit] <- substr(val[hit], 1L, pos[hit] - 1L)
    }
    data.table::set(dt, j = col, value = val)
  }

  if (nrow(dt) > 0L && anyDuplicated(dt$feature_id) > 0L) {
    rk <- if ("rank_final" %in% names(dt)) {
      suppressWarnings(as.numeric(dt$rank_final))
    } else {
      rep_len(1, nrow(dt))
    }
    sf <- -suppressWarnings(as.numeric(dt$score_final))
    dt <- dt[order(as.numeric(dt$feature_id), rk, sf, na.last = TRUE)]
  }

  dt[, `:=`(
    best_candidate_1 = gsub(" or .*", "", candidate_structure_tax_npc_01pat),
    best_candidate_2 = gsub(" or .*", "", candidate_structure_tax_npc_02sup),
    best_candidate_3 = gsub(" or .*", "", candidate_structure_tax_npc_03cla)
  )]

  out <- unique(
    dt[, list(
      feature_id = feature_id,
      mz = feature_mz,
      rt = feature_rt,
      smiles_2D = candidate_structure_smiles_no_stereo,
      inchikey_2D = candidate_structure_inchikey_connectivity_layer,
      molecular_formula = candidate_structure_molecular_formula,
      candidate_count_similarity_peaks_matched,
      candidate_score_similarity = candidate_score_similarity,
      score_initial = score_initial,
      score_biological = score_biological,
      score_chemical = score_chemical,
      score_final = score_final,
      best_candidate_organism = candidate_structure_organism_occurrence_closest,
      consensus_1 = feature_pred_tax_npc_01pat_val,
      consensus_2 = feature_pred_tax_npc_02sup_val,
      consensus_3 = feature_pred_tax_npc_03cla_val,
      consistency_1 = feature_pred_tax_npc_01pat_score,
      consistency_2 = feature_pred_tax_npc_02sup_score,
      consistency_3 = feature_pred_tax_npc_03cla_score,
      best_candidate_1 = best_candidate_1,
      best_candidate_2 = best_candidate_2,
      best_candidate_3 = best_candidate_3,
      mode = if ("mode" %in% names(dt)) mode else NA_character_
    )],
    by = c("feature_id", "mode")
  )

  for (col in names(out)) {
    data.table::set(out, j = col, value = y_as_na(out[[col]], ""))
  }

  na_smiles <- is.na(out$smiles_2D)
  out[na_smiles, best_candidate_1 := "notAnnotated"]
  out[na_smiles, best_candidate_2 := paste(best_candidate_1, "notAnnotated")]
  out[na_smiles, best_candidate_3 := paste(best_candidate_2, "notAnnotated")]

  has_smiles <- !na_smiles
  out[
    has_smiles & is.na(best_candidate_1),
    best_candidate_1 := "notClassified"
  ]
  out[
    has_smiles & is.na(best_candidate_2),
    best_candidate_2 := paste(best_candidate_1, "notClassified")
  ]
  out[
    has_smiles & is.na(best_candidate_3),
    best_candidate_3 := paste(best_candidate_2, "notClassified")
  ]

  out
}

#' Flag candidates below confidence score threshold
#'
#' @param df Data frame or `data.table` of candidates.
#' @param score Numeric confidence score threshold.
#'
#' @return A `data.table` with low-confidence annotations labeled
#'   `"notConfident"`.
#' @keywords internal
#' @noRd
make_confident <- function(df, score) {
  dt <- data.table::copy(data.table::as.data.table(df))
  sf <- suppressWarnings(as.numeric(dt$score_final))
  na_sf <- is.na(sf)
  low_sf <- !na_sf & sf < score
  dt[
    na_sf,
    `:=`(
      best_candidate_1 = NA_character_,
      best_candidate_2 = NA_character_,
      best_candidate_3 = NA_character_
    )
  ]
  dt[
    low_sf,
    `:=`(
      best_candidate_1 = "notConfident",
      best_candidate_2 = "notConfident notConfident",
      best_candidate_3 = "notConfident notConfident notConfident"
    )
  ]
  dt
}

#' Collapse minor chemical superclasses into "Other"
#'
#' @param dataframe Data frame or `data.table`.
#' @param value Name of numeric column to sum by category (default
#'   `"peak_area"`).
#'
#' @return A `data.table` with harmonized `"Other"` subcategories.
#' @keywords internal
#' @noRd
make_other <- function(dataframe, value = "peak_area") {
  dt <- data.table::copy(data.table::as.data.table(dataframe))
  if (nrow(dt) == 0L) {
    return(dt)
  }

  tmp <- data.table::copy(dt)
  tmp[,
    new := sum(get(value)),
    by = c("best_candidate_1", "best_candidate_2")
  ]
  data.table::setorder(tmp, -new, na.last = TRUE)
  tmp[, row_number := seq_len(.N), by = "best_candidate_1"]

  top_4 <- unique(tmp[
    row_number <= 4L,
    list(
      smiles_2D,
      inchikey_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3
    )
  ])

  last_dt <- dt[
    !top_4,
    on = c("best_candidate_1", "best_candidate_2")
  ]
  last <- unique(last_dt[, list(
    smiles_2D,
    inchikey_2D,
    best_candidate_1,
    best_candidate_2 = "Other",
    best_candidate_3
  )])

  new_map <- data.table::rbindlist(list(top_4, last))[,
    list(best_candidate_1, best_candidate_2, best_candidate_3)
  ]

  dt[, best_candidate_2 := NULL]
  unique(dt_left_join(dt, new_map))
}

#' Filter out unannotated, unclassified, and low-confidence rows
#'
#' @param dataframe Data frame or `data.table`.
#'
#' @return Filtered `data.table`.
#' @keywords internal
#' @noRd
no_other <- function(dataframe) {
  dt <- data.table::as.data.table(dataframe)
  dt[!grepl("not", best_candidate_1) & !is.na(best_candidate_1)]
}

#' Deduplicate stereoisomers to 2D structures
#'
#' @param df Data frame or `data.table` of structures.
#'
#' @return A `data.table` containing unique 2D structures per taxon/reference.
#' @keywords internal
#' @noRd
make_no_stereo <- function(df) {
  message("Keeping no stereo structures only")
  dt <- data.table::copy(data.table::as.data.table(df))
  non_struct_cols <- names(dt)[!grepl("structure", names(dt))]
  dedup_cols <- intersect(
    c(
      non_struct_cols,
      "structure_smiles_2D",
      "reference_wikidata",
      "reference_doi"
    ),
    names(dt)
  )
  out <- unique(dt, by = dedup_cols)
  if ("structure_smiles_2D" %in% names(out)) {
    out[, structure_smiles_2D := NULL]
  }
  out
}

#' Filter structures to chromatographiable mass and logP range
#'
#' @param df Data frame or `data.table` of structures.
#' @param mass_min Minimum exact mass (default `50`).
#' @param mass_max Maximum exact mass (default `1500`).
#' @param logp_min Minimum XLogP (default `-1`).
#' @param logp_max Maximum XLogP (default `6`).
#'
#' @return Filtered `data.table`.
#' @keywords internal
#' @noRd
make_chromatographiable <- function(
  df,
  mass_min = 50,
  mass_max = 1500,
  logp_min = -1,
  logp_max = 6
) {
  message("Keeping chromatographiable structures only")
  dt <- data.table::as.data.table(df)
  dt[
    structure_exact_mass >= mass_min &
      structure_exact_mass <= mass_max &
      structure_xlogp >= logp_min &
      structure_xlogp <= logp_max
  ]
}

#' Prepare peak comparison tables across similarity and taxonomy filters
#'
#' @param features_informed Path to informed features table (or `NULL`).
#' @param features_not_informed Path to non-informed features table (or `NULL`).
#' @param candidates_confident Confident candidates `data.table`.
#' @param min_similarity_prefilter Minimum similarity pre-filter
#'   (default `0.6`).
#' @param min_similarity_filter Minimum similarity filter (default `0.8`).
#' @param mode Ionization mode (`"pos"` or `"neg"`).
#' @param show_example Logical; use package example data if `TRUE`.
#' @param default_peak_area Default peak area assigned to features outside
#'   detected peaks (default `0.001`).
#' @param correct_acn Logical; apply acetonitrile gradient response correction
#'   to CAD peak areas via `predict_response()` (default `FALSE`).
#' @param acn_time Numeric; gradient time in minutes to reach 100% acetonitrile
#'   (default `30.5`).
#' @param dvol Numeric; dwell time offset in minutes (default `0.5`).
#'
#' @return A named list of 9 `data.table` objects representing successive
#'   filtering stages.
#' @keywords internal
#' @noRd
prepare_comparison <- function(
  features_informed = NULL,
  features_not_informed = NULL,
  candidates_confident,
  min_similarity_prefilter = 0.6,
  min_similarity_filter = 0.8,
  mode = "pos",
  show_example = FALSE,
  default_peak_area = 0.001,
  correct_acn = FALSE,
  acn_time = 30.5,
  dvol = 0.5
) {
  target_mode <- mode
  peaks_compared <- data.table::copy(data.table::as.data.table(
    load_features_informed(features_informed, show_example = show_example)
  ))
  peaks_compared[, mode := target_mode]

  peaks_outside <- data.table::copy(data.table::as.data.table(
    load_features_not_informed(
      features_not_informed,
      show_example = show_example
    )
  ))
  peaks_outside[, mode := target_mode]

  num_cols <- c(
    "peak_id",
    "peak_rt_min",
    "peak_rt_apex",
    "peak_rt_max",
    "peak_area",
    "feature_id",
    "feature_rt",
    "feature_mz",
    "feature_area",
    "comparison_score"
  )
  for (col in intersect(num_cols, names(peaks_outside))) {
    data.table::set(
      peaks_outside,
      j = col,
      value = suppressWarnings(as.numeric(peaks_outside[[col]]))
    )
  }

  peaks_all <- data.table::rbindlist(
    list(peaks_compared, peaks_outside),
    use.names = TRUE,
    fill = TRUE
  )

  message("joining compared peaks and candidates")
  cand_dt <- data.table::copy(data.table::as.data.table(candidates_confident))
  cand_dt[, feature_id := as.numeric(feature_id)]

  join_candidates <- function(df) {
    dt <- data.table::copy(data.table::as.data.table(df))
    dt[, feature_id := as.numeric(feature_id)]
    out <- dt_left_join(dt, cand_dt)
    out[, rt := as.numeric(rt)]
    out
  }

  add_id_cols <- function(dt) {
    dt <- data.table::copy(dt)
    dt[, `:=`(
      id = sample,
      integral = peak_area,
      intensity = feature_area
    )]
    dt
  }

  peaks_maj <- add_id_cols(join_candidates(peaks_compared))
  if (correct_acn && nrow(peaks_maj) > 0L) {
    peaks_maj[,
      acn := 5 +
        95 / (acn_time / pmin(acn_time, pmax(0.01, peak_rt_apex - dvol - 0.5)))
    ]
    peaks_maj[,
      peak_area := predict_response(acn = acn, peak_area = peak_area)
    ]
    peaks_maj[, `:=`(
      integral = peak_area,
      acn = NULL
    )]
  }

  peaks_min <- add_id_cols(join_candidates(peaks_outside))
  peaks_min[, `:=`(
    peak_rt_apex = rt,
    peak_area = default_peak_area
  )]
  for (col in c(
    "peak_id",
    "peak_rt_min",
    "peak_rt_max",
    "feature_rt",
    "feature_mz",
    "feature_area",
    "comparison_score",
    "integral",
    "intensity"
  )) {
    if (col %in% names(peaks_min)) {
      data.table::set(
        peaks_min,
        j = col,
        value = suppressWarnings(as.numeric(peaks_min[[col]]))
      )
    }
  }

  message("keeping peaks similarities above desired (pre-)threshold only")
  peaks_maj_precor <- peaks_maj[
    !is.na(comparison_score) & comparison_score >= min_similarity_prefilter
  ]

  peaks_min_precor <- data.table::rbindlist(
    list(
      add_id_cols(peaks_maj[!peaks_maj_precor, on = names(peaks_maj)]),
      peaks_min
    ),
    use.names = TRUE,
    fill = TRUE
  )

  message("keeping multiple features only if none was reported in the species")
  peaks_maj_precor_taxo <- data.table::copy(peaks_maj_precor)
  if (nrow(peaks_maj_precor_taxo) > 0L) {
    bco <- peaks_maj_precor_taxo$best_candidate_organism
    sp <- peaks_maj_precor_taxo$species
    valid_pair <- !is.na(bco) & !is.na(sp)
    taxo_val <- integer(nrow(peaks_maj_precor_taxo))
    if (any(valid_pair)) {
      taxo_val[valid_pair] <- as.integer(mapply(
        grepl,
        pattern = bco[valid_pair],
        x = sp[valid_pair],
        MoreArgs = list(fixed = TRUE),
        USE.NAMES = FALSE
      ))
    }
    taxo_2_val <- as.integer(valid_pair & bco == sp)

    peaks_maj_precor_taxo[, `:=`(
      taxo = as.numeric(taxo_val),
      taxo_2 = as.numeric(taxo_2_val)
    )]
    peaks_maj_precor_taxo[,
      `:=`(
        sum = sum(taxo),
        sum_2 = sum(taxo_2)
      ),
      by = c("sample", "peak_id")
    ]
    peaks_maj_precor_taxo[,
      keep := ifelse(
        taxo_2 == 1,
        "Y",
        ifelse(
          taxo == 1,
          ifelse(sum != sum_2 & sum_2 == 0, "Y", "N"),
          ifelse(sum == 0, "Y", "N")
        )
      )
    ]
    peaks_maj_precor_taxo <- peaks_maj_precor_taxo[keep == "Y"]
  } else {
    peaks_maj_precor_taxo[, `:=`(
      taxo = numeric(0),
      taxo_2 = numeric(0),
      sum = numeric(0),
      sum_2 = numeric(0),
      keep = character(0)
    )]
  }

  peaks_min_precor_taxo <- data.table::rbindlist(
    list(
      add_id_cols(peaks_maj[!peaks_maj_precor_taxo, on = names(peaks_maj)]),
      peaks_min
    ),
    use.names = TRUE,
    fill = TRUE
  )

  message(
    "keeping peaks similarities with score above ",
    min_similarity_filter
  )
  peaks_maj_precor_taxo_cor <- peaks_maj_precor_taxo[
    !is.na(comparison_score) & comparison_score >= min_similarity_filter
  ]

  peaks_min_precor_taxo_cor <- data.table::rbindlist(
    list(
      add_id_cols(peaks_maj[!peaks_maj_precor_taxo_cor, on = names(peaks_maj)]),
      peaks_min
    ),
    use.names = TRUE,
    fill = TRUE
  )

  list(
    peaks_all = peaks_all,
    peaks_maj = peaks_maj,
    peaks_min = peaks_min,
    peaks_maj_precor = peaks_maj_precor,
    peaks_min_precor = peaks_min_precor,
    peaks_maj_precor_taxo = peaks_maj_precor_taxo,
    peaks_min_precor_taxo = peaks_min_precor_taxo,
    peaks_maj_precor_taxo_cor = peaks_maj_precor_taxo_cor,
    peaks_min_precor_taxo_cor = peaks_min_precor_taxo_cor
  )
}
