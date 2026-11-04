library(tinytest)
library(cascade)

get_peaks <- utils::getFromNamespace("get_peaks", "cascade")

t_pts <- seq(1, 50, by = 1)
sig_1 <- 5 + 100 * exp(-((t_pts - 15)^2) / 8) + 80 * exp(-((t_pts - 35)^2) / 10)
sig_2 <- 2 + 50 * exp(-((t_pts - 25)^2) / 6)

mat_1col <- matrix(
  sig_1,
  ncol = 1L,
  dimnames = list(as.character(t_pts), "254")
)
mat_2col <- cbind("254" = sig_1, "280" = sig_2)
rownames(mat_2col) <- as.character(t_pts)

# 1. Matrix input (coerced to list) + missing lambdas on 1-column matrix +
# warning on unnamed list
expect_warning(
  res_mat <- get_peaks(mat_1col, fit = "gaussian", smooth_type = "none"),
  pattern = "Sample names not found"
)
expect_inherits(res_mat, "peak_list")
expect_identical(attr(res_mat, "lambdas"), "254")
expect_identical(attr(res_mat, "fit"), "gaussian")
expect_identical(attr(res_mat, "time.units"), "min")
expect_identical(nrow(res_mat[["1"]][["254"]]), 2L)

# 2. Error when lambdas is missing on multi-column matrix
expect_error(
  get_peaks(list(s1 = mat_2col)),
  pattern = "Wavelengths \\(`lambdas`\\) must be provided"
)

# 3. Error when specified wavelength does not exist
expect_error(
  get_peaks(list(s1 = mat_2col), lambdas = "999"),
  pattern = "could not be found"
)

# 4. Numeric lambdas coerced to character + collapse = TRUE
res_col <- get_peaks(
  list(s1 = mat_2col),
  lambdas = c(254, 280),
  fit = "raw",
  collapse = TRUE,
  smooth_type = "none"
)
expect_inherits(res_col$s1, "data.frame")
expect_identical(unique(res_col$s1$lambda), c("254", "280"))
expect_identical(nrow(res_col$s1), 3L)

# 5. Time unit scaling: "min" (1x), "s" (60x), "ms" (60000x)
res_min <- get_peaks(
  list(s1 = mat_1col),
  lambdas = "254",
  fit = "egh",
  time.units = "min",
  smooth_type = "none"
)
res_sec <- get_peaks(
  list(s1 = mat_1col),
  lambdas = "254",
  fit = "egh",
  time.units = "s",
  smooth_type = "none"
)
res_ms <- get_peaks(
  list(s1 = mat_1col),
  lambdas = "254",
  fit = "egh",
  time.units = "ms",
  smooth_type = "none"
)
expect_equal(
  res_sec$s1[["254"]]$area,
  res_min$s1[["254"]]$area * 60,
  tolerance = 1e-10
)
expect_equal(
  res_ms$s1[["254"]]$area,
  res_min$s1[["254"]]$area * 60000,
  tolerance = 1e-8
)
expect_equal(
  res_sec$s1[["254"]]$FWHM,
  res_min$s1[["254"]]$FWHM * 60,
  tolerance = 1e-10
)
expect_identical(res_sec$s1[["254"]]$rt, res_min$s1[["254"]]$rt)

# 6. All smoothing types in find_peaks: gaussian, box, savgol (even & odd
# window), mva, tmva, none
for (st in c("gaussian", "box", "savgol", "mva", "tmva", "none")) {
  r_st <- get_peaks(
    list(s1 = mat_1col),
    lambdas = "254",
    fit = "gaussian",
    smooth_type = st,
    smooth_window = 6
  )
  expect_true(nrow(r_st$s1[["254"]]) >= 1L)
}
# Savgol with odd smooth_window
r_sg_odd <- get_peaks(
  list(s1 = mat_1col),
  lambdas = "254",
  fit = "gaussian",
  smooth_type = "savgol",
  smooth_window = 5
)
expect_true(nrow(r_sg_odd$s1[["254"]]) >= 1L)

# 7. Flat signal -> no peaks detected -> returns 0-row data.frame with expected
# columns
mat_flat <- matrix(
  rep(1, 20L),
  ncol = 1L,
  dimnames = list(as.character(1:20), "254")
)
res_flat <- get_peaks(list(s1 = mat_flat), lambdas = "254", fit = "egh")
expect_identical(nrow(res_flat$s1[["254"]]), 0L)
expect_true(all(
  c("sample", "lambda", "rt", "start", "end", "area") %in%
    names(res_flat$s1[["254"]])
))

# 8. sd.max filtering removes wide peaks when small, keeps all when NULL
res_sd_strict <- get_peaks(
  list(s1 = mat_1col),
  lambdas = "254",
  fit = "gaussian",
  sd.max = 0.5,
  smooth_type = "none"
)
expect_identical(nrow(res_sd_strict$s1[["254"]]), 0L)

res_sd_null <- get_peaks(
  list(s1 = mat_1col),
  lambdas = "254",
  fit = "gaussian",
  sd.max = NULL,
  smooth_type = "none"
)
expect_identical(nrow(res_sd_null$s1[["254"]]), 2L)
