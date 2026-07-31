library(tinytest)

utils::data(
  cascade_annotations,
  cascade_features,
  cascade_features_informed,
  cascade_features_not_informed,
  cascade_chromatograms_positive,
  cascade_chromatograms_negative,
  cascade_ms_data,
  package = "cascade"
)

expect_true(is.data.frame(cascade_annotations))
expect_true(is.data.frame(cascade_features))
expect_true(is.data.frame(cascade_features_informed))
expect_true(is.data.frame(cascade_features_not_informed))
expect_true(is.data.frame(cascade_chromatograms_positive))
expect_true(is.data.frame(cascade_chromatograms_negative))
expect_true(inherits(cascade_ms_data, "MSnExp"))

expect_equal(
  names(cascade_chromatograms_positive),
  c("chromatogram", "rtime", "intensity")
)
expect_equal(
  names(cascade_chromatograms_negative),
  c("chromatogram", "rtime", "intensity")
)
expect_equal(
  unname(lapply(
    load_chromatograms(show_example = TRUE, example_polarity = "pos"),
    names
  )),
  list(
    c("rtime", "intensity"),
    c("rtime", "intensity"),
    c("rtime", "intensity")
  )
)
