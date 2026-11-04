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

# --- 2. prepare_tima_annotations with unnamed annotations vector ---
local({
  ann_tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(ann_tmp), add = TRUE)
  writeLines(
    c(
      "feature_id\tfeature_mz\tfeature_rt\tcandidate_structure_smiles_no_stereo\tcandidate_structure_inchikey_connectivity_layer\tcandidate_structure_molecular_formula\tcandidate_count_similarity_peaks_matched\tcandidate_score_similarity\tscore_initial\tscore_biological\tscore_chemical\tscore_final\tcandidate_structure_organism_occurrence_closest\tfeature_pred_tax_npc_01pat_val\tfeature_pred_tax_npc_02sup_val\tfeature_pred_tax_npc_03cla_val\tfeature_pred_tax_npc_01pat_score\tfeature_pred_tax_npc_02sup_score\tfeature_pred_tax_npc_03cla_score\tcandidate_structure_tax_npc_01pat\tcandidate_structure_tax_npc_02sup\tcandidate_structure_tax_npc_03cla\tfeature_spectrum_peaks\tcandidate_library",
      "1\t300.1\t5.0\tC=C\tIK1\tC2H4\t5\t0.9\t0.9\t0.8\t0.8\t0.85\tSwertia\tTerpenoids\tMono\tIridoids\t0.9\t0.9\t0.9\tTerpenoids\tMonoterpenoids\tSecoiridoids\t6\tgnps"
    ),
    ann_tmp
  )
  # Unnamed character vector triggers nms <- as.character(seq_along(annotations))
  tima_unnamed <- prepare_tima_annotations(
    annotations = ann_tmp,
    show_example = FALSE
  )
  expect_identical(names(tima_unnamed), "1")
  expect_identical(nrow(tima_unnamed[["1"]]), 1L)
})

# --- 3. generate_tables with file_negative, empty lotus match, and export = FALSE ---
local({
  wf_dir <- tempfile("cascade_wf_")
  dir.create(
    file.path(wf_dir, "data/source/libraries"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  on.exit(unlink(wf_dir, recursive = TRUE), add = TRUE)

  withr::with_dir(wf_dir, {
    # Create a tiny lotus.csv.gz with non-matching InChIKey to exercise nrow(lotus_sub) == 0L branch
    lotus_gz <- gzfile("data/source/libraries/lotus.csv.gz", open = "wb")
    writeLines(
      c(
        "structure_inchikey,organism_taxonomy_01domain,organism_taxonomy_02kingdom,organism_taxonomy_03phylum,organism_taxonomy_04class,organism_taxonomy_05order,organism_taxonomy_06family,organism_taxonomy_07tribe,organism_taxonomy_08genus,organism_taxonomy_09species,organism_taxonomy_10varietas,reference_doi",
        "NONMATCHINGIK-UHFFFAOYSA-N,Eukaryota,Plantae,Streptophyta,Magnoliopsida,Gentianales,Gentianaceae,Gentianeae,Swertia,Swertia chirayita,NA,10.1000/test"
      ),
      lotus_gz
    )
    close(lotus_gz)

    pos_tsv <- file.path(wf_dir, "pos.tsv")
    neg_tsv <- file.path(wf_dir, "neg.tsv")
    ann_tsv <- file.path(wf_dir, "ann.tsv")

    writeLines(
      c(
        "sample\tpeak_id\tpeak_rt_min\tpeak_rt_apex\tpeak_rt_max\tpeak_area\tfeature_id\tfeature_rt\tfeature_mz\tfeature_area\tcomparison_score",
        "s1\t1\t1.0\t1.5\t2.0\t100\t1\t1.5\t300.1\t1000\t0.95",
        "s1\t2\t3.0\t3.5\t4.0\t200\t2\t3.5\t400.2\t2000\t0.92"
      ),
      pos_tsv
    )
    writeLines(
      c(
        "sample\tpeak_id\tpeak_rt_min\tpeak_rt_apex\tpeak_rt_max\tpeak_area\tfeature_id\tfeature_rt\tfeature_mz\tfeature_area\tcomparison_score",
        "s1\t1\t1.0\t1.5\t2.0\t100\t1\t1.5\t298.1\t800\t0.90"
      ),
      neg_tsv
    )
    writeLines(
      c(
        "feature_id\tfeature_mz\tfeature_rt\tcandidate_structure_smiles_no_stereo\tcandidate_structure_inchikey_connectivity_layer\tcandidate_structure_molecular_formula\tcandidate_count_similarity_peaks_matched\tcandidate_score_similarity\tscore_initial\tscore_biological\tscore_chemical\tscore_final\tcandidate_structure_organism_occurrence_closest\tfeature_pred_tax_npc_01pat_val\tfeature_pred_tax_npc_02sup_val\tfeature_pred_tax_npc_03cla_val\tfeature_pred_tax_npc_01pat_score\tfeature_pred_tax_npc_02sup_score\tfeature_pred_tax_npc_03cla_score\tcandidate_structure_tax_npc_01pat\tcandidate_structure_tax_npc_02sup\tcandidate_structure_tax_npc_03cla",
        "1\t300.1\t1.5\tC=C\tIK1\tC2H4\t5\t0.9\t0.9\t0.8\t0.8\t0.85\tSwertia\tTerpenoids\tMono\tIridoids\t0.9\t0.9\t0.9\tTerpenoids\tMonoterpenoids\tSecoiridoids",
        "2\t400.2\t3.5\tC=O\tIK2\tCH2O\t6\t0.95\t0.95\t0.9\t0.9\t1.20\tSwertia\tAlkaloids\tIndole\tSimple\t0.9\t0.9\t0.9\tAlkaloids\tIndole\tSimple"
      ),
      ann_tsv
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
    expect_inherits(gt_res$csv_table, "data.table")
    expect_identical(nrow(gt_res$csv_table), 1L)
    expect_identical(gt_res$csv_table$`Annotation Score`, ">1.00")
    expect_false(dir.exists("data/processed"))
  })
})

# --- 4. generate_ids offline with comparison, show_title = FALSE, format_tables = TRUE ---
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
            "structure,structureLabel,structure_id,structureSmiles,taxaLabels,taxa,referencesLabels,references_ids,references",
            paste0("Q1,Comp1,IK1,C=C,", sp_prefix, " sp1,QT1,Ref1,10.1/1,QR1"),
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
        "structure_inchikey,structure_smiles_2D,structure_exact_mass,structure_xlogp,structure_taxonomy_npclassifier_01pathway,structure_taxonomy_npclassifier_02superclass,structure_taxonomy_npclassifier_03class",
        "IK1,C=C,300.0,1.5,Terpenoids,Monoterpenoids,Secoiridoids",
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
  })
})
