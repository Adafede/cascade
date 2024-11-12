library(testthat)

test_that(desc = "Bitter is better", code = {
  message("\nBitter is better.\n")
  succeed()
})
test_that(desc = "Check chromatogram alignment", code = {
  message("\n")
  check_chromatograms_alignment(show_example = TRUE)
  succeed()
})
test_that(desc = "Check peaks integration", code = {
  message("\n")
  check_peaks_integration(show_example = TRUE)
  succeed()
})
test_that(desc = "Process compare peaks", code = {
  message("\n")
  process_compare_peaks(show_example = TRUE)
  succeed()
})
test_that(desc = "Generate pseudochromatograms", code = {
  message("\n")
  generate_pseudochromatograms(show_example = TRUE)
  succeed()
})
test_that(desc = "Generate IDs", code = {
  message("\n")
  generate_ids()
  succeed()
})
test_that(desc = "Generate Tables", code = {
  message("\n")
  generate_tables()
  succeed()
})
