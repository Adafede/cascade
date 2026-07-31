library(tinytest)

message("Bitter is better.")
expect_true(TRUE)

message("")
expect_true({
  check_chromatograms_alignment(show_example = TRUE)
  TRUE
})

message("")
expect_true({
  check_peaks_integration(show_example = TRUE)
  TRUE
})

message("")
expect_true({
  process_compare_peaks(show_example = TRUE)
  TRUE
})

message("")
expect_true({
  generate_pseudochromatograms(show_example = TRUE)
  TRUE
})

message("")
expect_true({
  generate_tables(show_example = TRUE)
  TRUE
})

message("")
expect_true({
  generate_ids()
  TRUE
})
