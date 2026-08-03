library(tinytest)

# Run all tests in a temporary directory so that file-writing functions
# (e.g. generate_tables, process_compare_peaks) and the tima log file
# ("tima.log") do not pollute inst/tinytest/ or the package installation dir.
# withr::with_dir restores the working directory via on.exit(), even on error.
withr::with_dir(tempdir(), {
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
})
