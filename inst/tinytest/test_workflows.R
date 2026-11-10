library(tinytest)
library(cascade)

add_chromato_line <- utils::getFromNamespace("add_chromato_line", "cascade")
check_chromatograms <- utils::getFromNamespace("check_chromatograms", "cascade")

# --- 1. add_chromato_line & check_chromatograms (neg polarity + validation) ---
df_tr <- data.frame(rtime = c(1, 2, 4), intensity = c(10, 50, 25))
p_base <- plotly::plot_ly()
p_neg <- add_chromato_line(
  p_base,
  chromato = df_tr,
  shift = 0,
  normalize_time = FALSE,
  name = "neg_trace",
  color = "#000000",
  polarity = "neg"
)
expect_inherits(p_neg, "plotly")

expect_error(
  check_chromatograms(
    chromatograms = "invalid_key",
    chromatograms_list = list()
  ),
  pattern = "chromatograms must be in"
)
expect_error(
  check_chromatograms(
    chromatograms = "cad_pos",
    chromatograms_list = list(),
    type = "raw"
  ),
  pattern = "type must be one of"
)

p_all_neg <- check_chromatograms(
  chromatograms = c("bpi_neg", "cad_neg", "pda_neg"),
  chromatograms_list = list(
    chromatogram_bpi_neg_baselined = df_tr,
    chromatogram_cad_neg_baselined = df_tr,
    chromatogram_pda_neg_baselined = df_tr
  ),
  normalize_time = TRUE,
  type = "baselined"
)
expect_inherits(p_all_neg, "plotly")

p_raw_time <- check_chromatograms(
  chromatograms = "bpi_pos",
  chromatograms_list = list(chromatogram_bpi_pos_baselined = df_tr),
  normalize_time = FALSE,
  type = "baselined"
)
expect_inherits(p_raw_time, "plotly")

# check_chromatograms_alignment with auto shift, neg-only, no-BPI, and file paths
local({
  ca_ns <- asNamespace("cascade")
  orig_load_chrom <- get("load_chromatograms", envir = ca_ns)
  unlockBinding("load_chromatograms", ca_ns)
  on.exit(
    {
      assign("load_chromatograms", orig_load_chrom, envir = ca_ns)
      lockBinding("load_chromatograms", ca_ns)
    },
    add = TRUE
  )
  t_g <- seq(0.5, 2.0, by = 0.02)
  mock_chroms <- list(
    bpi = data.frame(
      rtime = t_g,
      intensity = 10 + 100 * exp(-((t_g - 1.2)^2) / 0.01)
    ),
    pda = data.frame(
      rtime = t_g - 0.02,
      intensity = 5 + 80 * exp(-((t_g - 1.2)^2) / 0.01)
    ),
    cad = data.frame(
      rtime = t_g - 0.04,
      intensity = 2 + 90 * exp(-((t_g - 1.2)^2) / 0.01)
    )
  )
  assign(
    "load_chromatograms",
    function(...) mock_chroms,
    envir = ca_ns
  )

  # 1. Pos + neg with cad_shift = "auto", pda_shift = "auto" and file paths
  p_auto_pos <- check_chromatograms_alignment(
    file_positive = "pos.mzML",
    file_negative = "neg.mzML",
    show_example = FALSE,
    improve_signal = FALSE,
    baseline_method = "als",
    cad_shift = "auto",
    pda_shift = "auto"
  )
  expect_inherits(p_auto_pos, "plotly")

  # 2. Negative-only chromatograms with cad_shift = "auto", pda_shift = "auto"
  p_auto_neg <- check_chromatograms_alignment(
    file_negative = "neg.mzML",
    chromatograms = c("bpi_neg", "cad_neg", "pda_neg"),
    show_example = FALSE,
    improve_signal = FALSE,
    baseline_method = "als",
    cad_shift = "auto",
    pda_shift = "auto"
  )
  expect_inherits(p_auto_neg, "plotly")

  # 3. No BPI selected -> ref_bpi is NULL -> fallback shift 0
  p_auto_nobpi <- check_chromatograms_alignment(
    file_positive = "pos.mzML",
    chromatograms = "cad_pos",
    show_example = FALSE,
    improve_signal = FALSE,
    baseline_method = "als",
    cad_shift = "auto",
    pda_shift = "auto"
  )
  expect_inherits(p_auto_nobpi, "plotly")
})

# Package startup message (.onAttach)
on_attach <- utils::getFromNamespace(".onAttach", "cascade")
expect_silent(suppressPackageStartupMessages(on_attach("lib", "cascade")))

# --- 2. prepare_tima_annotations with unnamed annotations vector ---
local({
  ann_tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(ann_tmp), add = TRUE)
  ann_df <- data.table::data.table(
    feature_id = 1L,
    feature_mz = 300.1,
    feature_rt = 5.0,
    candidate_structure_smiles_no_stereo = "C=C",
    candidate_structure_inchikey_connectivity_layer = "IK1",
    candidate_structure_molecular_formula = "C2H4",
    candidate_count_similarity_peaks_matched = 5L,
    candidate_score_similarity = 0.9,
    score_initial = 0.9,
    score_biological = 0.8,
    score_chemical = 0.8,
    score_final = 0.85,
    candidate_structure_organism_occurrence_closest = "Swertia",
    feature_pred_tax_npc_01pat_val = "Terpenoids",
    feature_pred_tax_npc_02sup_val = "Mono",
    feature_pred_tax_npc_03cla_val = "Iridoids",
    feature_pred_tax_npc_01pat_score = 0.9,
    feature_pred_tax_npc_02sup_score = 0.9,
    feature_pred_tax_npc_03cla_score = 0.9,
    candidate_structure_tax_npc_01pat = "Terpenoids",
    candidate_structure_tax_npc_02sup = "Monoterpenoids",
    candidate_structure_tax_npc_03cla = "Secoiridoids",
    feature_spectrum_peaks = 6L,
    candidate_library = "gnps"
  )
  data.table::fwrite(ann_df, ann_tmp, sep = "\t")
  # Unnamed character vector triggers sequential integer names
  tima_unnamed <- prepare_tima_annotations(
    annotations = ann_tmp,
    show_example = FALSE
  )
  expect_identical(names(tima_unnamed), "1")
  expect_identical(nrow(tima_unnamed[["1"]]), 1L)
})

# --- 3. generate_tables with file_negative, empty lotus match, export=FALSE ---
local({
  wf_dir <- tempfile("cascade_wf_")
  dir.create(
    file.path(wf_dir, "data/source/libraries"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  on.exit(unlink(wf_dir, recursive = TRUE), add = TRUE)

  withr::with_dir(wf_dir, {
    # Tiny lotus.csv.gz with non-matching InChIKey to exercise empty match
    lotus_gz <- gzfile("data/source/libraries/lotus.csv.gz", open = "wb")
    writeLines(
      c(
        paste(
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
          "reference_doi",
          sep = ","
        ),
        paste(
          "NONMATCHINGIK-UHFFFAOYSA-N",
          "Eukaryota",
          "Plantae",
          "Streptophyta",
          "Magnoliopsida",
          "Gentianales",
          "Gentianaceae",
          "Gentianeae",
          "Swertia",
          "Swertia chirayita",
          "NA",
          "10.1000/test",
          sep = ","
        )
      ),
      lotus_gz
    )
    close(lotus_gz)

    pos_tsv <- file.path(wf_dir, "pos.tsv")
    neg_tsv <- file.path(wf_dir, "neg.tsv")
    ann_tsv <- file.path(wf_dir, "ann.tsv")

    data.table::fwrite(
      data.table::data.table(
        sample = c("s1", "s1", "s1", "s1"),
        peak_id = c(1L, 2L, 3L, 4L),
        peak_rt_min = c(1.0, 3.49, 5.0, 3.52),
        peak_rt_apex = c(1.5, 3.50, NA_real_, 3.52),
        peak_rt_max = c(2.0, 3.51, 6.0, 3.53),
        peak_area = c(100, 200, 50, 150),
        feature_id = c(1L, 2L, 1L, 1L),
        feature_rt = c(1.5, 3.50, 1.5, 3.52),
        feature_mz = c(300.1, 400.2, 300.1, 300.1),
        feature_area = c(1000, 2000, 500, 900),
        comparison_score = c(0.95, 0.92, 0.91, 0.90)
      ),
      pos_tsv,
      sep = "\t"
    )
    data.table::fwrite(
      data.table::data.table(
        sample = "s1",
        peak_id = 1L,
        peak_rt_min = 1.0,
        peak_rt_apex = 1.5,
        peak_rt_max = 2.0,
        peak_area = 100,
        feature_id = 1L,
        feature_rt = 1.5,
        feature_mz = 298.1,
        feature_area = 800,
        comparison_score = 0.90
      ),
      neg_tsv,
      sep = "\t"
    )
    data.table::fwrite(
      data.table::data.table(
        feature_id = c(1L, 2L),
        feature_mz = c(300.1, 400.2),
        feature_rt = c(1.5, 3.5),
        candidate_structure_smiles_no_stereo = c("C=C", "C=O"),
        candidate_structure_inchikey_connectivity_layer = c("IK1", "IK2"),
        candidate_structure_molecular_formula = c("C2H4", "CH2O"),
        candidate_count_similarity_peaks_matched = c(5L, 6L),
        candidate_score_similarity = c(0.9, 0.95),
        score_initial = c(0.9, 0.95),
        score_biological = c(0.8, 0.9),
        score_chemical = c(0.8, 0.9),
        score_final = c(0.85, 1.20),
        candidate_structure_organism_occurrence_closest = c(
          "Swertia",
          "Swertia"
        ),
        feature_pred_tax_npc_01pat_val = c("Terpenoids", "Alkaloids"),
        feature_pred_tax_npc_02sup_val = c("Mono", "Indole"),
        feature_pred_tax_npc_03cla_val = c("Iridoids", "Simple"),
        feature_pred_tax_npc_01pat_score = c(0.9, 0.9),
        feature_pred_tax_npc_02sup_score = c(0.9, 0.9),
        feature_pred_tax_npc_03cla_score = c(0.9, 0.9),
        candidate_structure_tax_npc_01pat = c("Terpenoids", "Alkaloids"),
        candidate_structure_tax_npc_02sup = c("Monoterpenoids", "Indole"),
        candidate_structure_tax_npc_03cla = c("Secoiridoids", "Simple")
      ),
      ann_tsv,
      sep = "\t"
    )

    gt_res <- generate_tables(
      annotations = ann_tsv,
      file_positive = pos_tsv,
      file_negative = neg_tsv,
      show_example = FALSE,
      export_csv = FALSE,
      export_html = FALSE
    )
    expect_inherits(gt_res$pretty_table, "gt_tbl")
    expect_equal(
      gt_res$csv_table[["Peak RT [min]"]],
      c(1.5, 3.5),
      tolerance = 1e-12
    )
    expect_identical(as.integer(gt_res$csv_table[["Feature ID"]]), c(1L, 2L))
    expect_inherits(gt_res$csv_table, "data.table")
    expect_identical(nrow(gt_res$csv_table), 2L)
    expect_identical(gt_res$csv_table$`Annotation Score`, c("0.85", ">1.00"))
    expect_false(dir.exists("data/processed"))
  })
})

# --- 4. generate_ids offline with comparison and format_tables = TRUE ---
local({
  ids_dir <- tempfile("cascade_ids_")
  dir.create(
    file.path(ids_dir, "data/source/libraries"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  on.exit(unlink(ids_dir, recursive = TRUE), add = TRUE)

  curl_ns <- asNamespace("curl")
  orig_dl <- get("curl_download", envir = curl_ns)
  unlockBinding("curl_download", curl_ns)
  on.exit(
    {
      assign("curl_download", orig_dl, envir = curl_ns)
      lockBinding("curl_download", curl_ns)
    },
    add = TRUE
  )

  state <- new.env(parent = emptyenv())
  state$q_call <- 0L
  assign(
    "curl_download",
    function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
      state$q_call <- state$q_call + 1L
      if (state$q_call <= 2L) {
        # taxon_name_to_qid calls for Swertia and Gentiana
        qid <- if (state$q_call == 1L) "Q158572" else "Q158500"
        writeLines(
          c("taxon", paste0("http://www.wikidata.org/entity/", qid)),
          destfile
        )
      } else {
        # wiki_progress calls
        sp_prefix <- if (state$q_call == 3L) "Swertia" else "Gentiana"
        writeLines(
          c(
            paste(
              "structure",
              "structureLabel",
              "structure_id",
              "structureSmiles",
              "taxaLabels",
              "taxa",
              "referencesLabels",
              "references_ids",
              "references",
              sep = ","
            ),
            paste0("Q1,Comp1,IK1,C=C,", sp_prefix, " sp1,QT1,Ref1,10.1/1,QR1"),
            paste0(
              "Q1b,Comp1b,IK1b,C=C,",
              sp_prefix,
              " sp1,QT1,Ref1,10.1/1,QR1"
            ),
            paste0("Q2,Comp2,IK2,C=O,", sp_prefix, " sp2,QT2,Ref2,10.1/2,QR2")
          ),
          destfile
        )
      }
      invisible(destfile)
    },
    envir = curl_ns
  )

  withr::with_dir(ids_dir, {
    lotus_gz <- gzfile("data/source/libraries/lotus.csv.gz", open = "wb")
    writeLines(
      c(
        paste(
          "structure_inchikey",
          "structure_smiles_2D",
          "structure_exact_mass",
          "structure_xlogp",
          "structure_taxonomy_npclassifier_01pathway",
          "structure_taxonomy_npclassifier_02superclass",
          "structure_taxonomy_npclassifier_03class",
          sep = ","
        ),
        "IK1,C=C,300.0,1.5,Terpenoids,Monoterpenoids,Secoiridoids",
        "IK1b,C=C,300.0,1.5,Terpenoids,Monoterpenoids,Secoiridoids",
        "IK2,C=O,350.0,2.0,Alkaloids,Indole alkaloids,Simple indole"
      ),
      lotus_gz
    )
    close(lotus_gz)

    res_ids <- generate_ids(
      taxa = c("Swertia", "Gentiana"),
      comparison = c("Swertia", "Gentiana"),
      no_stereo = FALSE,
      filter_ms_conditions = FALSE,
      show_title = FALSE,
      format_tables = TRUE
    )
    expect_true(all(
      c("plots", "treemaps", "sunbursts", "tables") %in% names(res_ids)
    ))
    expect_true("special" %in% names(res_ids$treemaps))
    expect_inherits(res_ids$tables$Swertia, "gt_tbl")
    expect_inherits(res_ids$tables$Gentiana, "gt_tbl")
    expect_identical(nrow(res_ids$tables$Swertia[["_data"]]), 3L)
  })
})

# --- 5. check_peaks_integration(chromatogram = "original") and process_compare_peaks(0 peaks) ---
local({
  ca_ns <- asNamespace("cascade")
  orig_load_chrom <- get("load_chromatograms", envir = ca_ns)
  orig_load_ms <- get("load_ms_data", envir = ca_ns)
  unlockBinding("load_chromatograms", ca_ns)
  unlockBinding("load_ms_data", ca_ns)
  on.exit(
    {
      assign("load_chromatograms", orig_load_chrom, envir = ca_ns)
      assign("load_ms_data", orig_load_ms, envir = ca_ns)
      lockBinding("load_chromatograms", ca_ns)
      lockBinding("load_ms_data", ca_ns)
    },
    add = TRUE
  )
  assign(
    "load_chromatograms",
    function(...) orig_load_chrom(show_example = TRUE),
    envir = ca_ns
  )
  assign(
    "load_ms_data",
    function(...) {
      list(
        rt = c(60, 90, 120),
        mz = list(300.1, 300.1, 300.1),
        intensity = list(100, 500, 100)
      )
    },
    envir = ca_ns
  )

  feat_tmp <- tempfile(fileext = ".csv")
  exp_dir <- tempfile("cascade_exp_")
  on.exit(unlink(c(feat_tmp, exp_dir), recursive = TRUE), add = TRUE)
  feat_df <- data.frame(
    id = 1L,
    rt = 99.0,
    mz = 300.1,
    area = 100000,
    status = "DETECTED",
    rt_min = 98.8,
    rt_max = 99.2,
    mz_min = 300.09,
    mz_max = 300.11,
    int_min = 100,
    int_max = 50000,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(feat_df)[5:11] <- c(
    "datafile:s1.mzML:feature_state",
    "datafile:s1.mzML:rt_range:min",
    "datafile:s1.mzML:rt_range:max",
    "datafile:s1.mzML:mz_range:min",
    "datafile:s1.mzML:mz_range:max",
    "datafile:s1.mzML:intensity_range:min",
    "datafile:s1.mzML:intensity_range:max"
  )
  data.table::fwrite(feat_df, feat_tmp)

  p_orig <- check_peaks_integration(
    file = "s1.mzML",
    features = feat_tmp,
    chromatogram = "original",
    show_example = FALSE
  )
  expect_inherits(p_orig, "plotly")

  res_empty_comp <- process_compare_peaks(
    file = "s1.mzML",
    features = feat_tmp,
    export_dir = exp_dir,
    show_example = FALSE
  )
  expect_identical(nrow(res_empty_comp$features_informed), 0L)
  expect_identical(nrow(res_empty_comp$features_not_informed), 1L)
})
