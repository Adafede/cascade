library(tinytest)
library(cascade)

compare_peaks_fn <- utils::getFromNamespace("compare_peaks", "cascade")
get_peaks_fn <- utils::getFromNamespace("get_peaks", "cascade")

# Regression test 1: compare_peaks() must return 0 (not crash with NaN
# tolerance)
# when the detector peak has <= 1 retention time point.
peaks_prelist_single <- list(
  list_chromato_peaks = list(
    list(rtime = 0.5, intensity = 1.0),
    list(rtime = numeric(0L), intensity = numeric(0L))
  )
)
list_ms_peaks_multi <- list(
  list(list(rtime = c(0.0, 0.5, 1.0), intensity = c(0.2, 1.0, 0.3))),
  list(list(rtime = c(0.0, 0.5, 1.0), intensity = c(0.2, 1.0, 0.3)))
)
expect_identical(
  compare_peaks_fn(1L, list_ms_peaks_multi, peaks_prelist_single),
  list(0)
)
expect_identical(
  compare_peaks_fn(2L, list_ms_peaks_multi, peaks_prelist_single),
  list(0)
)

# Regression test 2: get_peaks(lambdas = "max") on a multi-column matrix must
# select the column with the global maximum rather than indexing by linear
# offset.
rt_grid <- seq(1, 10, by = 1)
col_254 <- c(0, 1, 5, 1, 0, 0, 0, 0, 0, 0)
col_280 <- c(0, 0, 0, 0, 0, 2, 15, 2, 0, 0)
mat_multi <- cbind("254" = col_254, "280" = col_280)
rownames(mat_multi) <- as.character(rt_grid)

res_max <- get_peaks_fn(
  list(sample1 = mat_multi),
  lambdas = "max",
  fit = "raw",
  smooth_type = "none"
)
expect_inherits(res_max, "peak_list")
expect_identical(res_max$sample1$max$rt, 6)
expect_identical(res_max$sample1$max$height, 2)
expect_identical(res_max$sample1$max$start, 4)
expect_identical(res_max$sample1$max$end, 8)

# Regression test 3: get_peaks() on a flat chromatogram with 0 detected peaks
# must return a 0-row data.frame instead of coercing noPeaksMat to a character
# matrix via cbind() and failing with 'non-numeric argument to binary operator'.
mat_flat <- matrix(
  rep(1, 20L),
  ncol = 1L,
  dimnames = list(as.character(1:20), "254")
)
res_flat <- get_peaks_fn(list(s1 = mat_flat), lambdas = "254", fit = "egh")
expect_inherits(res_flat$s1[["254"]], "data.frame")
expect_identical(nrow(res_flat$s1[["254"]]), 0L)

# Regression test 4: get_peaks() when sd.max filters out all detected peaks
# must return a 0-row data.frame instead of failing in cbind() with
# 'arguments imply differing number of rows: 1, 0'.
res_sd_filtered <- get_peaks_fn(
  list(s1 = mat_multi),
  lambdas = "280",
  fit = "gaussian",
  sd.max = 0.01,
  smooth_type = "none"
)
expect_inherits(res_sd_filtered$s1[["280"]], "data.frame")
expect_identical(nrow(res_sd_filtered$s1[["280"]]), 0L)

# Regression test 5: filter_fft() must preserve the DC component when
# round(n * components) < 2L and must not reverse the index range when
# components >= 0.5 or n <= 2L.
filter_fft_fn <- utils::getFromNamespace("filter_fft", "cascade")
expect_equal(
  filter_fft_fn(rep(10, 50L), components = 0.01),
  rep(10, 50L),
  tolerance = 1e-12
)
expect_equal(
  filter_fft_fn(c(2, 4, 6, 8), components = 1.0),
  c(2, 4, 6, 8),
  tolerance = 1e-12
)
expect_equal(
  filter_fft_fn(c(3, 7), components = 0.01),
  c(3, 7),
  tolerance = 1e-12
)
expect_equal(
  filter_fft_fn(c(1, 5, -1, -5), components = 0.01),
  c(0.5, 2.5, -0.5, -2.5),
  tolerance = 1e-12
)
