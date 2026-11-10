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

# 9. Exact recovery of tailed EGH parameters (sigma, tau, FWHM, area, R^2 = 1)
t_egh <- seq(1, 100, by = 1)
true_sd <- 3.5
true_tau <- 2.5
true_h <- 100
denom_e <- 2 * true_sd^2 + true_tau * (t_egh - 45)
y_egh_test <- ifelse(denom_e > 0, true_h * exp(-(t_egh - 45)^2 / denom_e), 0)
mat_egh_test <- matrix(
  y_egh_test,
  ncol = 1L,
  dimnames = list(as.character(t_egh), "254")
)
res_egh_exact <- get_peaks(
  list(s1 = mat_egh_test),
  lambdas = "254",
  fit = "egh",
  smooth_type = "none"
)$s1[["254"]]
expect_equal(res_egh_exact$sd[1L], true_sd, tolerance = 1e-3)
expect_equal(res_egh_exact$tau[1L], true_tau, tolerance = 1e-3)
expect_equal(
  res_egh_exact$FWHM[1L],
  sqrt(8 * log(2) * true_sd^2 + (log(2) * true_tau)^2),
  tolerance = 1e-3
)
expect_equal(res_egh_exact$r.squared[1L], 1.0, tolerance = 1e-6)
expect_equal(
  res_egh_exact$area[1L],
  sum((y_egh_test[-1L] + y_egh_test[-length(y_egh_test)]) / 2),
  tolerance = 1e-3
)
expect_true(res_min$s1[["254"]]$r.squared[1L] > 0.99)

# 10. Zero-height peak (start_height == 0) yields ssy == 0 -> r.squared == 0
mat_zero_h <- matrix(
  c(-10, -8, -5, 0, -5, -8, -10),
  ncol = 1L,
  dimnames = list(as.character(1:7), "254")
)
res_zero_h <- get_peaks(
  list(s1 = mat_zero_h),
  lambdas = "254",
  fit = "gaussian",
  smooth_type = "none",
  amp_thresh = -20,
  noise_threshold = -1
)$s1[["254"]]
expect_identical(nrow(res_zero_h), 1L)
expect_equal(res_zero_h$r.squared[1L], 0, tolerance = 1e-12)

# 11. High-baseline peak where no point in [start, end] drops below half_h
# exercises fallback t_left / t_right and max.iter = 0L
mat_high_base <- matrix(
  c(10, 8, 9, 10, 9, 8, 10),
  ncol = 1L,
  dimnames = list(as.character(1:7), "254")
)
res_high_base <- get_peaks(
  list(s1 = mat_high_base),
  lambdas = "254",
  fit = "egh",
  max.iter = 0L,
  smooth_type = "none"
)$s1[["254"]]
expect_identical(nrow(res_high_base), 1L)
expect_true(res_high_base$area[1L] > 0)

# 12. lambdas = "max" selects wavelength column with highest maximum intensity
res_max_lam <- get_peaks(
  list(s1 = mat_2col),
  lambdas = "max",
  fit = "gaussian",
  smooth_type = "none"
)$s1[["max"]]
expect_identical(nrow(res_max_lam), 2L)
expect_equal(res_max_lam$rt, c(15, 35), tolerance = 0.1)
