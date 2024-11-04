library(testthat)

## need to do all in one because of outputs needed in the same temp dir
## use fixtures instead in the future
test_that(desc = "Test functions", code = {
  testthat::expect_no_error("Yayy" != "Yaay")
  message("Bitter is better.")
  check_chromatograms_alignment(show_example = TRUE)
  check_peaks_integration(show_example = TRUE)
  process_compare_peaks(show_example = TRUE)
  process_plot_pseudochromatograms(show_example = TRUE)
  succeed()
})
