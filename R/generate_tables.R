#' Generate interactive HTML and TSV summary tables of annotated peaks
#'
#' Integrates peak-informed features across positive and/or negative modes with
#' TIMA structural annotations and LOTUS taxonomic/literature metadata, applies
#' confidence-ordered 1D non-maximum suppression across retention-time windows,
#' and produces both an interactive HTML table (with 2D chemical structure
#' depictions via CDK Depict) and a tab-separated summary table.
#'
#' @export
#'
#' @param annotations Path to TIMA annotations `.tsv` file (optional when
#'   `show_example = TRUE`).
#' @param file_negative Path to negative-mode informed features `.tsv` file
#'   (optional).
#' @param file_positive Path to positive-mode informed features `.tsv` file
#'   (optional when `show_example = TRUE`).
#' @param min_confidence Minimum final annotation confidence score. Default is
#'   `0.4`.
#' @param show_example Logical; use bundled example data. Default is `FALSE`.
#' @param export_csv Logical; export `.tsv` summary table. Default is `TRUE`.
#' @param export_html Logical; export `.html` interactive `gt` table. Default
#'   is `TRUE`.
#' @param export_dir Directory where output files are written. Default is
#'   `"data/processed"`.
#' @param export_name Base filename (without extension) for exported tables.
#'   Default is `"cascade_table"`.
#'
#' @return A named list with two elements:
#'   - `pretty_table`: A `gt_tbl` interactive table object.
#'   - `csv_table`: A `data.table` containing the formatted summary rows.
#'
#' @concept EDAM:operation_3803 Natural product identification
#' @concept EDAM:data_0962 Small molecule report
#' @concept EDAM:data_1712 Chemical structure image
#' @concept EDAM:format_2331 HTML
#' @concept EDAM:format_3475 TSV
#' @concept EDAM:topic_3172 Metabolomics
#' @concept EDAM:topic_0154 Small molecules
#'
#' @examples
#' \donttest{
#' old_wd <- setwd(tempdir())
#' generate_tables(show_example = TRUE, export_csv = FALSE, export_html = FALSE)
#' setwd(old_wd)
#' }
generate_tables <- function(
  annotations = NULL,
  file_negative = NULL,
  file_positive = NULL,
  min_confidence = 0.4,
  show_example = FALSE,
  export_csv = TRUE,
  export_html = TRUE,
  export_dir = "data/processed",
  export_name = "cascade_table"
) {
  if (!show_example && is.null(file_positive) && is.null(file_negative)) {
    stop(
      "At least one of 'file_positive' or 'file_negative' must be provided.",
      call. = FALSE
    )
  }

  message("loading informed features")
  tables_peaks_pos <- NULL
  tables_peaks_neg <- NULL

  if (!is.null(file_positive) || show_example) {
    message("... positive mode")
    tables_peaks_pos <- data.table::as.data.table(
      if (show_example) {
        load_features_informed(show_example = TRUE)
      } else {
        load_features_informed(file_positive)
      }
    )
  }
  if (!is.null(file_negative)) {
    message("... negative mode")
    neg_dt <- data.table::as.data.table(
      load_features_informed(file_negative, show_example = show_example)
    )
    tables_peaks_neg <- unique(neg_dt)
  }

  message("loading annotations")
  ann_mode <- if (is.null(file_positive) && !show_example) "neg" else "pos"
  annotation_table <- keep_best_candidates(
    load_annotations(annotations, show_example = show_example, mode = ann_mode)
  )
  annotation_table[, feature_id := as.numeric(feature_id)]
  annotation_table <- make_confident(annotation_table, score = min_confidence)

  lotus_path <- "data/source/libraries/lotus.csv.gz"
  message("Getting last LOTUS version")
  get_zenodo_file(
    doi = "10.5281/zenodo.5794106",
    pattern = "frozen_metadata.csv.gz",
    path = lotus_path
  )

  message("Loading LOTUS classified structures")
  lotus_cols <- c(
    "structure_inchikey",
    "organism_taxonomy_01domain",
    "organism_taxonomy_02kingdom",
    "organism_taxonomy_03phylum",
    "organism_taxonomy_04class",
    "organism_taxonomy_05order",
    "organism_taxonomy_06family",
    "organism_taxonomy_07tribe",
    "organism_taxonomy_08genus",
    "organism_taxonomy_09species",
    "organism_taxonomy_10varietas",
    "reference_doi"
  )
  lotus <- unique(fread_gz(
    file = lotus_path,
    select = lotus_cols,
    colClasses = "character"
  ))

  peak_parts <- Filter(
    Negate(is.null),
    list(tables_peaks_pos, tables_peaks_neg)
  )
  table_peaks <- if (length(peak_parts) > 1L) {
    data.table::rbindlist(peak_parts, use.names = TRUE, fill = TRUE)
  } else {
    data.table::copy(peak_parts[[1L]])
  }

  data.table::setkey(table_peaks, peak_rt_min, peak_rt_max)

  table_joined <- dt_inner_join(table_peaks, annotation_table)
  table_processed <- table_joined[, list(
    peak_rt_min = peak_rt_min,
    peak_rt_max = peak_rt_max,
    feature_id = feature_id,
    feature_rt = feature_rt,
    feature_mz = feature_mz,
    smiles_no_stereo = smiles_2D,
    best_candidate_1 = best_candidate_1,
    best_candidate_2 = best_candidate_2,
    best_candidate_3 = best_candidate_3,
    score_final = score_final,
    best_candidate_organism = best_candidate_organism,
    inchikey_connectivity_layer = inchikey_2D
  )]

  data.table::setkey(table_processed, peak_rt_min, peak_rt_max)

  message("joining within given rt tolerance")
  overlapped <- data.table::foverlaps(table_peaks, table_processed)
  data.table::setorder(overlapped, -score_final, na.last = TRUE)
  overlapped <- unique(overlapped, by = c("peak_rt_apex", "smiles_no_stereo"))
  overlapped[,
    `:=`(
      newrt = mean(peak_rt_apex),
      newarea = mean(peak_area)
    ),
    by = "peak_rt_min"
  ]
  overlapped <- unique(overlapped, by = c("newrt", "newarea"))
  if (nrow(overlapped) > 0L) {
    sf_num <- suppressWarnings(as.numeric(overlapped$score_final))
    ord_score <- order(
      -sf_num,
      -overlapped$newarea,
      overlapped$newrt,
      na.last = TRUE
    )
    overlapped <- overlapped[ord_score]
    keep_rows <- logical(nrow(overlapped))
    kept_rts <- numeric(0L)
    for (i in seq_len(nrow(overlapped))) {
      rt_i <- overlapped$newrt[i]
      if (is.na(rt_i)) {
        next
      }
      if (length(kept_rts) == 0L || all(abs(rt_i - kept_rts) > 0.05)) {
        keep_rows[i] <- TRUE
        kept_rts <- c(kept_rts, rt_i)
      }
    }
    overlapped <- overlapped[keep_rows]
    data.table::setorder(overlapped, newrt, na.last = TRUE)
  }
  overlapped[, peak_id := seq_len(.N)]

  total_area <- sum(overlapped$newarea, na.rm = TRUE)
  pct_area <- if (is.finite(total_area) && total_area > 0) {
    round(100 * overlapped$newarea / total_area, digits = 2)
  } else {
    rep(0, nrow(overlapped))
  }

  table_medium <- overlapped[, list(
    peak_id = peak_id,
    peak_rt = round(newrt, digits = 2),
    peak_area = pct_area,
    feature_id = feature_id,
    feature_rt = round(feature_rt, digits = 2),
    feature_mz = round(feature_mz, digits = 4),
    smiles_no_stereo = smiles_no_stereo,
    best_candidate_1 = best_candidate_1,
    best_candidate_2 = best_candidate_2,
    best_candidate_3 = best_candidate_3,
    score_final = round(as.numeric(score_final), digits = 2),
    best_candidate_organism = best_candidate_organism,
    inchikey_connectivity_layer = inchikey_connectivity_layer
  )]

  lotus_sub <- data.table::copy(lotus)
  lotus_sub[, inchikey_connectivity_layer := sub("-.*", "", structure_inchikey)]
  target_iks <- unique(table_medium$inchikey_connectivity_layer)
  lotus_sub <- lotus_sub[inchikey_connectivity_layer %in% target_iks]

  tax_cols <- lotus_cols[2:11]
  if (nrow(lotus_sub) > 0L) {
    lotus_long <- unique(data.table::melt(
      lotus_sub,
      id.vars = c("inchikey_connectivity_layer", "reference_doi"),
      measure.vars = tax_cols,
      variable.name = "rank",
      value.name = "best_candidate_organism",
      na.rm = FALSE
    )[, list(
      inchikey_connectivity_layer,
      reference = reference_doi,
      best_candidate_organism
    )])
  } else {
    lotus_long <- data.table::data.table(
      inchikey_connectivity_layer = character(0),
      reference = character(0),
      best_candidate_organism = character(0)
    )
  }

  tf_join <- dt_left_join(table_medium, lotus_long)
  tf_join[,
    reference := {
      refs <- unique(
        reference[!is.na(reference) & nzchar(reference) & reference != "NA"]
      )
      if (length(refs) == 0L) NA_character_ else paste(refs, collapse = "; ")
    },
    by = c("peak_id", "inchikey_connectivity_layer")
  ]
  tf_join[,
    Structure := encode_smiles(smiles_no_stereo, empty = "", reserved = TRUE)
  ]
  tf_join[,
    annotation_score := ifelse(
      score_final > 1,
      ">1.00",
      as.character(score_final)
    )
  ]

  table_final <- unique(tf_join[, list(
    `Peak ID` = peak_id,
    `Peak RT [min]` = peak_rt,
    `Peak Area [%]` = peak_area,
    `Feature ID` = feature_id,
    `Feature RT [min]` = feature_rt,
    `Feature m/z` = feature_mz,
    Structure = Structure,
    `InChIKey no stereo` = inchikey_connectivity_layer,
    `Annotation Score` = annotation_score,
    `Chemical Pathway` = best_candidate_1,
    `Chemical Superclass` = best_candidate_2,
    `Chemical Class` = best_candidate_3,
    `Closest Organism with Structure Reported` = best_candidate_organism,
    `Reference(s)` = reference
  )])

  table_pretty <- table_final |>
    gt::gt(auto_align = TRUE) |>
    gt::text_transform(
      locations = gt::cells_body(columns = "Structure"),
      fn = molinfo
    ) |>
    gt::sub_missing(
      columns = gt::everything(),
      rows = gt::everything(),
      missing_text = ""
    ) |>
    gt::opt_interactive(use_filters = TRUE)

  if (export_csv || export_html) {
    message("checking export directory")
    check_export_dir(export_dir)
    message("exporting")
  }

  if (export_csv) {
    message("...CSV")
    data.table::fwrite(
      table_final,
      file = file.path(export_dir, paste0(export_name, ".tsv")),
      sep = "\t"
    )
  }
  if (export_html) {
    message("...HTML")
    gt::gtsave(
      table_pretty,
      filename = file.path(export_dir, paste0(export_name, ".html"))
    )
  }

  list("pretty_table" = table_pretty, "csv_table" = table_final)
}
