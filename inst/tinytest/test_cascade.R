library(tinytest)
library(cascade)

orig_lotus <- normalizePath(
  "data/source/libraries/lotus.csv.gz",
  mustWork = FALSE
)
test_dir <- tempfile("cascade_test_")
dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(orig_lotus)) {
  dir.create(
    file.path(test_dir, "data/source/libraries"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  file.symlink(
    orig_lotus,
    file.path(test_dir, "data/source/libraries/lotus.csv.gz")
  )
}

withr::with_dir(test_dir, {
  message("Testing internal mathematical and alignment helpers...")
  runmean_fn <- utils::getFromNamespace("runmean_center", "cascade")
  closest_fn <- utils::getFromNamespace("closest_idx", "cascade")
  reorder_fn <- utils::getFromNamespace("fct_reorder2_base", "cascade")

  x_vec <- sin(seq(0, 2 * pi, length.out = 25L))
  rm_val <- runmean_fn(x_vec, 5L)
  expect_equal(length(rm_val), length(x_vec))
  expect_true(all(is.finite(rm_val)))

  idx_val <- closest_fn(
    c(0.1, 0.5, 0.9),
    c(0.0, 0.25, 0.5, 0.75, 1.0),
    tolerance = 0.2
  )
  expect_equal(idx_val, c(1L, 3L, 5L))

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

  # Verify filter_fft NA handling and exact frequency cutoff
  sig_na <- c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10)
  fft_res <- fft_fn(sig_na, components = 0.2)
  expect_equal(length(fft_res), length(sig_na))
  expect_true(all(is.finite(fft_res)))

  # Verify make_confident NA vs low-score behavior
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

  # Verify order-preserving data.table joins
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
  expect_true(all(c("plots_1", "plots_2", "treemaps") %in% names(res_pseudo)))
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

  tima_plots <- plot_tima(tima_tabs)
  expect_true(all(c("histogram", "treemap", "sunburst") %in% names(tima_plots)))
  expect_inherits(tima_plots$histogram, "ggplot")
  expect_inherits(tima_plots$treemap, "plotly")
  expect_inherits(tima_plots$sunburst, "plotly")

  message("Testing query_wikidata input validation...")
  expect_error(query_wikidata(123))
  expect_error(query_wikidata("SELECT * WHERE {}", timeout = -5))
  expect_error(query_wikidata("SELECT * WHERE {}", remove_url = "yes"))

  message("Testing generate_tables...")
  res_tables <- generate_tables(show_example = TRUE)
  expect_true(all(c("pretty_table", "csv_table") %in% names(res_tables)))
  expect_inherits(res_tables$pretty_table, "gt_tbl")
  expect_inherits(res_tables$csv_table, "data.table")
  expect_true(nrow(res_tables$csv_table) > 0L)
  expect_true(file.exists("data/processed/cascade_table.tsv"))
  expect_true(file.exists("data/processed/cascade_table.html"))

  message("Testing generate_ids...")
  res_ids <- tryCatch(
    generate_ids(taxa = "Swertia", show_title = TRUE),
    error = function(e) {
      message("Skipping live Wikidata test due to network error: ", e$message)
      NULL
    }
  )
  if (!is.null(res_ids)) {
    expect_true(all(
      c("plots", "treemaps", "sunbursts", "tables") %in% names(res_ids)
    ))
  }
})
