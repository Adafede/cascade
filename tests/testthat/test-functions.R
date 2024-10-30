library(testthat)

## need to do all in one because of outputs needed in the same temp dir
## use fixtures instead in the future
test_that(desc = "Test functions", code = {
  testthat::expect_no_error("Yayy" != "Yaay")
  message("Bitter is better.")
  succeed()
})
