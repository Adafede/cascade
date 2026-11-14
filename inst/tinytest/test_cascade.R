library(tinytest)
library(cascade)

local({
  test_dir <- tempfile("cascade_test_")
  dir.create(
    file.path(test_dir, "data/source/libraries"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  withr::local_envvar(
    c(R_USER_CACHE_DIR = file.path(test_dir, "user_cache"))
  )
  withr::local_dir(test_dir)

  # Create deterministic offline LOTUS fixture covering generate_tables() and
  # generate_ids() columns so tests never download 86 MB from Zenodo.
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
        "KOHAOVSIVQWVAQ-UHFFFAOYSA-N,C=C,300.1,1.5",
        "Terpenoids,Monoterpenoids,Secoiridoids",
        "Eukaryota,Plantae,Streptophyta,Magnoliopsida",
        "Gentianales,Gentianaceae,Gentianeae,Swertia",
        "Swertia chirayita,NA,10.1021/test.lotus1",
        sep = ","
      ),
      paste(
        "JZRWCGZRTZMZEH-UHFFFAOYSA-N,CC=O,320.2,2.0",
        "Terpenoids,Monoterpenoids,Secoiridoids",
        "Eukaryota,Plantae,Streptophyta,Magnoliopsida",
        "Gentianales,Gentianaceae,Gentianeae,Swertia",
        "Swertia perennis,NA,10.1021/test.lotus2",
        sep = ","
      ),
      paste(
        "VFRROHXSMXFLSN-UHFFFAOYSA-N,O=C1OCC2C1C(O)OC=C2,374.1,0.5",
        "Terpenoids,Monoterpenoids,Secoiridoids",
        "Eukaryota,Plantae,Streptophyta,Magnoliopsida",
        "Gentianales,Gentianaceae,Gentianeae,Swertia",
        "Swertia chirayita,NA,10.1021/test.lotus3",
        sep = ","
      )
    ),
    lotus_gz
  )
  close(lotus_gz)

  message("Testing internal mathematical and alignment helpers...")
  runmean_fn <- utils::getFromNamespace("runmean_center", "cascade")
  reorder_fn <- utils::getFromNamespace("fct_reorder2_base", "cascade")

  x_vec <- sin(seq(0, 2 * pi, length.out = 25L))
  rm_val <- runmean_fn(x_vec, 5L)
  expect_equal(length(rm_val), length(x_vec))
  expect_true(all(is.finite(rm_val)))

  f_ord <- reorder_fn(
    c("a", "b", "a", "b"),
    c(1, 2, 3, 4),
    c(10, 5, 20, 15),
    desc = TRUE
  )
  expect_equal(levels(f_ord), c("a", "b"))

  pred_fn <- utils::getFromNamespace("predict_response", "cascade")
  pal_list <- utils::getFromNamespace("cascade_palettes", "cascade")
  format_gt_fn <- utils::getFromNamespace("format_gt", "cascade")
  fft_fn <- utils::getFromNamespace("filter_fft", "cascade")
  conf_fn <- utils::getFromNamespace("make_confident", "cascade")
  lj_fn <- utils::getFromNamespace("dt_left_join", "cascade")
  ij_fn <- utils::getFromNamespace("dt_inner_join", "cascade")
  fj_fn <- utils::getFromNamespace("dt_full_join", "cascade")
  expect_equal(length(pal_list), 15L)
  resp_val <- pred_fn(acn = c(20, 50, 80), peak_area = c(0, 100, 1000))
  expect_equal(length(resp_val), 3L)
  expect_equal(resp_val[1L], 0)
  expect_true(all(resp_val[2:3] > 0))

  sig_na <- c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10)
  fft_res <- fft_fn(sig_na, components = 0.2)
  expect_equal(length(fft_res), length(sig_na))
  expect_true(all(is.finite(fft_res)))

  conf_in <- data.table::data.table(
    score_final = c(0.8, 0.2, NA_real_),
    best_candidate_1 = c("Alkaloids", "Terpenoids", "Shikimates"),
    best_candidate_2 = c("Indole", "Monoterpenoids", "Flavonoids"),
    best_candidate_3 = c("Simple indole", "Iridoids", "Flavones")
  )
  conf_out <- conf_fn(conf_in, score = 0.5)
  expect_equal(
    conf_out$best_candidate_1,
    c("Alkaloids", "notConfident", NA_character_)
  )
  expect_equal(
    conf_out$best_candidate_2,
    c("Indole", "notConfident notConfident", NA_character_)
  )
  expect_equal(
    conf_out$best_candidate_3,
    c("Simple indole", "notConfident notConfident notConfident", NA_character_)
  )

  dt_a <- data.table::data.table(k = c("b", "a", "c"), v1 = c(2L, 1L, 3L))
  dt_b <- data.table::data.table(k = c("a", "b"), v2 = c(10L, 20L))
  lj_res <- lj_fn(dt_a, dt_b, by = "k")
  expect_equal(lj_res$k, c("b", "a", "c"))
  expect_equal(names(lj_res), c("k", "v1", "v2"))
  expect_equal(lj_res$v2, c(20L, 10L, NA_integer_))
  ij_res <- ij_fn(dt_a, dt_b, by = "k")
  expect_equal(ij_res$k, c("a", "b"))
  fj_res <- fj_fn(dt_a, dt_b, by = "k")
  expect_equal(nrow(fj_res), 3L)

  sample_gt_df <- data.frame(
    structureLabel = "Swertiamarin",
    structure = "Q415869",
    structureImage = "O=C1OCC2C1C(O)OC=C2",
    taxaLabels = "Swertia chirayita",
    taxa = "Q10758325",
    referencesLabels = "Test paper",
    references = "Q123456",
    stringsAsFactors = FALSE
  )
  expect_inherits(format_gt_fn(sample_gt_df, title = "Swertia"), "gt_tbl")

  message("Testing check_chromatograms_alignment...")
  p_align <- check_chromatograms_alignment(show_example = TRUE)
  expect_inherits(p_align, "plotly")
  expect_inherits(p_align, "htmlwidget")

  p_align_raw <- check_chromatograms_alignment(
    show_example = TRUE,
    improve_signal = FALSE,
    normalize_time = TRUE,
    type = "improved"
  )
  expect_inherits(p_align_raw, "plotly")

  message("Testing check_peaks_integration...")
  p_peaks <- check_peaks_integration(show_example = TRUE)
  expect_inherits(p_peaks, "plotly")

  message("Testing process_compare_peaks...")
  res_comp <- process_compare_peaks(show_example = TRUE)
  expect_true(is.list(res_comp))
  expect_true(all(
    c("features_informed", "features_not_informed") %in% names(res_comp)
  ))
  expect_inherits(res_comp$features_informed, "data.table")
  expect_inherits(res_comp$features_not_informed, "data.table")
  expect_equal(nrow(res_comp$features_informed), 4132L)
  expect_equal(nrow(res_comp$features_not_informed), 650L)
  expect_true(
    file.exists(
      "data/interim/peaks/210619_AR_06_V_03_2_01_featuresInformed_cad.tsv"
    )
  )
  expect_true(
    file.exists(
      "data/interim/peaks/210619_AR_06_V_03_2_01_featuresNotInformed_cad.tsv"
    )
  )

  message("Testing generate_pseudochromatograms...")
  res_pseudo <- generate_pseudochromatograms(
    show_example = TRUE,
    correct_acn = TRUE
  )
  expect_true(is.list(res_pseudo))
  expect_true(all(
    c("plots_1", "plots_2", "treemaps") %in% names(res_pseudo)
  ))
  expect_equal(length(res_pseudo$plots_1), 8L)
  expect_equal(length(res_pseudo$plots_2), 2L)
  expect_inherits(res_pseudo$plots_1$histograms_taxo_maj, "ggplot")
  expect_inherits(res_pseudo$plots_2$sunburst_signal_based_duo, "plotly")
  expect_inherits(res_pseudo$treemaps$peaks_maj, "plotly")

  message("Testing prepare_tima_annotations and plot_tima...")
  tima_tabs <- prepare_tima_annotations(show_example = TRUE)
  expect_true(is.list(tima_tabs))
  expect_inherits(tima_tabs[[1L]], "data.table")
  expect_true(nrow(tima_tabs[[1L]]) > 0L)

  tima_pred <- prepare_tima_annotations(
    show_example = TRUE,
    predicted_classes = TRUE
  )
  expect_inherits(tima_pred[[1L]], "data.table")
  expect_identical(names(tima_pred), "annotations_pred")

  tima_plots <- plot_tima(tima_tabs)
  expect_true(all(
    c("histogram", "treemap", "sunburst") %in% names(tima_plots)
  ))
  expect_inherits(tima_plots$histogram, "ggplot")
  expect_inherits(tima_plots$treemap, "plotly")
  expect_inherits(tima_plots$sunburst, "plotly")

  message("Testing query_wikidata input validation...")
  expect_error(query_wikidata(123))
  expect_error(query_wikidata("SELECT * WHERE {}", timeout = -5))
  expect_error(query_wikidata("SELECT * WHERE {}", remove_url = "yes"))

  message("Testing generate_tables offline with LOTUS fixture...")
  res_tables <- generate_tables(show_example = TRUE)
  expect_true(all(c("pretty_table", "csv_table") %in% names(res_tables)))
  expect_inherits(res_tables$pretty_table, "gt_tbl")
  expect_inherits(res_tables$csv_table, "data.table")
  expect_true(nrow(res_tables$csv_table) > 0L)
  expect_true(file.exists("data/processed/cascade_table.tsv"))
  expect_true(file.exists("data/processed/cascade_table.html"))

  message("Testing generate_ids offline with mocked SPARQL...")
  local({
    curl_ns <- asNamespace("curl")
    orig_dl <- get("curl_download", envir = curl_ns)
    unlockBinding("curl_download", curl_ns)
    on.exit(
      {
        if (bindingIsLocked("curl_download", curl_ns)) {
          unlockBinding("curl_download", curl_ns)
        }
        assign("curl_download", orig_dl, envir = curl_ns)
        lockBinding("curl_download", curl_ns)
      },
      add = TRUE
    )
    state <- new.env(parent = emptyenv())
    state$calls <- 0L
    make_mock_dl <- function(st) {
      force(st)
      function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
        st$calls <- st$calls + 1L
        if (st$calls == 1L) {
          writeLines(
            c("taxon", "http://www.wikidata.org/entity/Q158572"),
            destfile
          )
        } else {
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
              paste(
                "Q415869,Swertiamarin,VFRROHXSMXFLSN-UHFFFAOYSA-N",
                "O=C1OCC2C1C(O)OC=C2,Swertia chirayita,Q10758325",
                "Ref1,10.1021/test.lotus3,QR1",
                sep = ","
              )
            ),
            destfile
          )
        }
        invisible(destfile)
      }
    }
    assign("curl_download", make_mock_dl(state), envir = curl_ns)

    res_ids <- generate_ids(taxa = "Swertia", show_title = TRUE)
    expect_true(all(
      c("plots", "treemaps", "sunbursts", "tables") %in% names(res_ids)
    ))
    expect_inherits(res_ids$plots$Swertia, "ggplot")
    expect_inherits(res_ids$treemaps$Swertia, "plotly")
    expect_inherits(res_ids$sunbursts$Swertia, "plotly")
    expect_true("Swertiamarin" %in% res_ids$tables$Swertia$structureLabel)
  })
})
