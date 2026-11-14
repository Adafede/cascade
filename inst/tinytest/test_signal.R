library(tinytest)
library(cascade)

runmean_center <- utils::getFromNamespace("runmean_center", "cascade")
deriv_fn <- utils::getFromNamespace("deriv", "cascade")
middle_pts <- utils::getFromNamespace("middle_pts", "cascade")
second_der <- utils::getFromNamespace("second_der", "cascade")
filter_fft <- utils::getFromNamespace("filter_fft", "cascade")
signal_sharpening <- utils::getFromNamespace("signal_sharpening", "cascade")
improve_signal <- utils::getFromNamespace("improve_signal", "cascade")
improve_signals_progress <- utils::getFromNamespace(
  "improve_signals_progress",
  "cascade"
)
baseline_chromatogram <- utils::getFromNamespace(
  "baseline_chromatogram",
  "cascade"
)
normalize_chromatograms_list <- utils::getFromNamespace(
  "normalize_chromatograms_list",
  "cascade"
)
change_intensity_name <- utils::getFromNamespace(
  "change_intensity_name",
  "cascade"
)
extract_chromatogram <- utils::getFromNamespace(
  "extract_chromatogram",
  "cascade"
)
preprocess_chromatograms <- utils::getFromNamespace(
  "preprocess_chromatograms",
  "cascade"
)

# --- runmean_center ---
expect_identical(runmean_center(numeric(0L), 3L), numeric(0L))
expect_identical(runmean_center(c(1, 2, 3), 1L), c(1, 2, 3))
expect_identical(runmean_center(c(1, 2, 3), 0L), c(1, 2, 3))
expect_identical(runmean_center(c(1, 2, 3), -2L), c(1, 2, 3))

# k == 2 branch: y[1] = x[1], y[i] = (x[i-1] + x[i]) / 2
expect_identical(
  runmean_center(c(10, 20, 40, 80), 2L),
  c(10, 15, 30, 60)
)

# Odd k (k = 3): k1 = 1, k2 = 1
expect_identical(
  runmean_center(c(3, 6, 9, 12, 15), 3L),
  c(4.5, 6, 9, 12, 13.5)
)

# Even k > 2 (k = 4): k2 = 2, k1 = 1
expect_identical(
  runmean_center(c(2, 4, 6, 8, 10), 4L),
  c(4, 5, 7, 8, 9)
)

# k > n clamps k to n
expect_identical(
  runmean_center(c(2, 4, 6), 10L),
  runmean_center(c(2, 4, 6), 3L)
)
expect_identical(runmean_center(c(2, 4, 6), 10L), c(3, 4, 5))

# Property: affine linearity of centered running mean
set.seed(42L)
x_rand <- rnorm(30L)
expect_equal(
  runmean_center(3.5 * x_rand - 7, 5L),
  3.5 * runmean_center(x_rand, 5L) - 7,
  tolerance = 1e-12
)

# --- deriv, middle_pts, second_der ---
expect_identical(deriv_fn(c(1, 3, 4), c(2, 8, 6)), c(3, -2))
expect_identical(middle_pts(c(1, 5, 6)), c(3, 5.5))
# Second derivative of y = 3 * x^2 on uniform grid is identically 6
x_grid <- as.numeric(1:6)
expect_equal(
  second_der(x_grid, 3 * x_grid^2),
  rep(6, 4L),
  tolerance = 1e-12
)
# Non-uniform grid check for second_der
expect_equal(
  second_der(c(0, 1, 3), c(0, 1, 9)),
  (4 - 1) / (2 - 0.5),
  tolerance = 1e-12
)

# --- filter_fft ---
x_fft <- c(2, NA, 4, 6, 8, 10, 8, 6, 4, 2)
x_zero_na <- c(2, 0, 4, 6, 8, 10, 8, 6, 4, 2)
expect_identical(
  filter_fft(x_fft, components = 0.2),
  filter_fft(x_zero_na, components = 0.2)
)
# Independent reference calculation (reflection-padded 4th-order Butterworth)
x_pad <- c(x_zero_na, rev(x_zero_na))
n_pad <- length(x_pad)
cutoff_fft <- max(2, n_pad * 0.3)
freq_idx <- pmin(seq_len(n_pad) - 1L, n_pad - seq_len(n_pad) + 1L)
bw_w <- 1 / (1 + (freq_idx / cutoff_fft)^8)
ref_fft <- (Re(stats::fft(stats::fft(x_pad) * bw_w, inverse = TRUE)) / n_pad)[
  seq_along(x_zero_na)
]
expect_equal(
  filter_fft(x_zero_na, components = 0.3),
  ref_fft,
  tolerance = 1e-12
)

# --- signal_sharpening ---
runmean_symmetric <- utils::getFromNamespace("runmean_symmetric", "cascade")
estimate_rt_shift <- utils::getFromNamespace("estimate_rt_shift", "cascade")
expect_identical(runmean_symmetric(numeric(0L), 4L), numeric(0L))
expect_identical(runmean_symmetric(c(1, 2, 3), 1L), c(1, 2, 3))
expect_identical(
  runmean_symmetric(c(3, 6, 9, 12, 15), 3L),
  runmean_center(c(3, 6, 9, 12, 15), 3L)
)
expect_identical(signal_sharpening(1:4, 1:4), numeric(0L))

t_sharp <- seq(0, 10, length.out = 101L)
y_const <- rep(5, length(t_sharp))
# On constant signal, all derivatives are 0, so output == 5 + Baseline_adjust
expect_equal(
  signal_sharpening(t_sharp, y_const, Baseline_adjust = 2.5),
  rep(7.5, length(t_sharp) - 4L),
  tolerance = 1e-12
)
# With sigma = 0, sharpening equals first-pass symmetric smoothed signal on
# interior indices 3:(n - 2)
y_gauss <- exp(-((t_sharp - 5)^2) / 0.5)
s0 <- signal_sharpening(t_sharp, y_gauss, sigma = 0, Smoothing_width = 4L)
s_pos <- signal_sharpening(
  t_sharp,
  y_gauss,
  sigma = 0.05,
  k2 = 10,
  k4 = 100,
  Smoothing_width = 4L
)
expect_identical(length(s0), length(t_sharp) - 4L)
expect_equal(
  s0,
  runmean_symmetric(y_gauss, k = 4L)[3:(length(t_sharp) - 2L)],
  tolerance = 1e-12
)
expect_true(max(s_pos) > max(s0))
# Positive 4th-deriv term at Gaussian apex boosts height over 2nd-deriv only
s_no_d4 <- signal_sharpening(
  t_sharp,
  y_gauss,
  sigma = 0.05,
  k2 = 10,
  k4 = Inf,
  Smoothing_width = 4L
)
expect_true(s_pos[49L] > s_no_d4[49L])
# Symmetry preservation: symmetric Gaussian centered at midpoint index 51 (in
# 1..101 -> index 49 in 3..99) remains strictly symmetric
expect_equal(s_pos[49L - (1:20)], s_pos[49L + (1:20)], tolerance = 1e-12)

# --- estimate_rt_shift ---
t_shift_grid <- seq(0, 5, by = 0.01)
y_shift_ref <- exp(-((t_shift_grid - 2.5)^2) / 0.02)
df_shift_ref <- data.frame(rtime = t_shift_grid, intensity = y_shift_ref)
df_shift_tgt <- data.frame(
  rtime = t_shift_grid - 0.06,
  intensity = y_shift_ref * 0.5
)
expect_equal(
  estimate_rt_shift(
    df_shift_ref,
    df_shift_tgt,
    max_shift = 0.2,
    grid_step = 0.01
  ),
  0.06,
  tolerance = 1e-12
)
expect_identical(
  estimate_rt_shift(
    data.frame(rtime = c(1, 2), intensity = c(5, 5)),
    data.frame(rtime = c(1, 2), intensity = c(5, 5))
  ),
  0
)
expect_identical(
  estimate_rt_shift(
    data.frame(rtime = c(1, 2), intensity = c(1, 2)),
    data.frame(rtime = c(5, 6), intensity = c(1, 2))
  ),
  0
)

# --- improve_signal and improve_signals_progress ---
df_pos <- data.frame(
  rtime = seq(0, 2, length.out = 50L),
  intensity = seq(1, 10, length.out = 50L)
)
df_neg <- data.frame(
  rtime = seq(0, 2, length.out = 50L),
  intensity = seq(-5, 4, length.out = 50L)
)
df_zero <- data.frame(
  rtime = seq(0, 2, length.out = 50L),
  intensity = seq(0, 9, length.out = 50L)
)
imp_pos <- improve_signal(
  df_pos,
  fourier_components = 0.1,
  frequency = 1,
  resample = 1,
  time_min = 0,
  time_max = 1.5,
  intensity_floor = 0.5
)
imp_pos_floor2 <- improve_signal(
  df_pos,
  fourier_components = 0.1,
  frequency = 1,
  resample = 1,
  time_min = 0,
  time_max = 1.5,
  intensity_floor = 5.0
)
imp_zero <- improve_signal(
  df_zero,
  fourier_components = 0.1,
  frequency = 1,
  resample = 1,
  time_min = 0,
  time_max = 1.5,
  intensity_floor = 1.0
)
imp_neg <- improve_signal(
  df_neg,
  fourier_components = 0.1,
  frequency = 1,
  resample = 1,
  time_min = 0,
  time_max = 1.5,
  intensity_floor = 6.0
)
# Expected time grid: interior points seq(0, 1.5, by = 1/60)[3:89]
expected_time <- seq(0, 1.5, by = 1 / 60)[3:89]
expect_equal(imp_pos$rtime, expected_time, tolerance = 1e-12)
expect_equal(imp_neg$rtime, expected_time, tolerance = 1e-12)
# When min_intensity > 0, shift_amount is 0 regardless of intensity_floor
expect_equal(imp_pos$intensity, imp_pos_floor2$intensity, tolerance = 1e-12)
# When min_intensity == 0, shift_amount = 0 + 1.0 = 1.0, matching df_pos
expect_equal(imp_zero$intensity, imp_pos$intensity, tolerance = 1e-12)
# When min_intensity is -5, shift_amount is 11.0, matching df_pos shifted by 5
expect_equal(imp_neg$intensity, imp_pos$intensity + 5.0, tolerance = 1e-12)

prog_res <- improve_signals_progress(
  list(a = df_pos, b = df_neg),
  fourier_components = 0.1,
  frequency = 1,
  resample = 1,
  time_min = 0,
  time_max = 1.5
)
expect_identical(names(prog_res), c("a", "b"))
expect_equal(prog_res$a, imp_pos, tolerance = 1e-12)

# --- baseline_chromatogram ---
df_bl <- data.frame(
  rtime = seq(1, 20, by = 1),
  intensity = c(
    10,
    NA,
    10,
    12,
    30,
    80,
    30,
    12,
    10,
    10,
    10,
    10,
    15,
    40,
    15,
    10,
    10,
    10,
    10,
    10
  )
)
bl_out <- baseline_chromatogram(df_bl, method = "als")
expect_identical(bl_out$rtime, df_bl$rtime)
expect_identical(nrow(bl_out), nrow(df_bl))
expect_true(all(is.finite(bl_out$intensity)))

# --- normalize_chromatograms_list ---
raw_chrom <- data.frame(rtime = c(2, 4, 8), intensity = c(10, 50, 25))
n_int_only <- normalize_chromatograms_list(
  raw_chrom,
  shift = 2,
  normalize_intensity = TRUE,
  normalize_time = FALSE
)
expect_equal(
  n_int_only,
  data.frame(rtime = c(4, 6, 10), intensity = c(0.2, 1.0, 0.5)),
  tolerance = 1e-12
)

n_time_only <- normalize_chromatograms_list(
  raw_chrom,
  shift = 2,
  normalize_intensity = FALSE,
  normalize_time = TRUE
)
expect_equal(
  n_time_only,
  data.frame(
    rtime = c(0.4, 0.6, 1.0),
    intensity = c(10, 50, 25),
    time_2 = c(10, 10, 10)
  ),
  tolerance = 1e-12
)

n_neither <- normalize_chromatograms_list(
  raw_chrom,
  shift = -1,
  normalize_intensity = FALSE,
  normalize_time = FALSE
)
expect_equal(
  n_neither,
  data.frame(rtime = c(1, 3, 7), intensity = c(10, 50, 25)),
  tolerance = 1e-12
)

# --- change_intensity_name & extract_chromatogram ---
custom_df <- data.frame(ret_t = c(1.5, 2.5), sig = c(100, 200))
expect_identical(
  change_intensity_name(custom_df, name_rt = "ret_t", name_intensity = "sig"),
  data.frame(rtime = c(1.5, 2.5), intensity = c(100, 200))
)

hdrs <- c("bpi" = "BPI_H", "cad" = "CAD_H", "pda" = "PDA_H")
chrom_bundle <- list(
  bpi = data.frame(rtime = c(1, 2), intensity = c(10, 20)),
  cad = data.frame(rtime = c(1, 2), intensity = c(30, 40)),
  pda = data.frame(rtime = c(1, 2), intensity = c(50, 60))
)
expect_identical(
  extract_chromatogram(chrom_bundle, "cad", hdrs),
  data.frame(rtime = c(1, 2), intensity = c(30, 40))
)
expect_error(
  extract_chromatogram(chrom_bundle, "uv", hdrs),
  pattern = "should be one of"
)

# 2-detector bundle missing pda: CAD still resolves by name, PDA errors clearly
two_det <- chrom_bundle[c("bpi", "cad")]
expect_identical(
  extract_chromatogram(two_det, "cad", hdrs),
  data.frame(rtime = c(1, 2), intensity = c(30, 40))
)
expect_error(
  extract_chromatogram(two_det, "pda", hdrs),
  pattern = "not found"
)

# Interior NAs in runmean_center must not poison adjacent values
rm_na <- runmean_center(c(10, NA_real_, 30, 40, 50), k = 3L)
expect_false(anyNA(rm_na))
expect_equal(rm_na[2L], 20, tolerance = 1e-12)
rm_all_na <- runmean_center(c(NA_real_, NA_real_, NA_real_), k = 3L)
expect_true(all(is.na(rm_all_na)))

# improve_signal with solvent delay (min(rtime) > time_min) and short grid
delay_df <- data.frame(
  rtime = seq(1.5, 4.0, length.out = 200L),
  intensity = 10 +
    100 *
      exp(
        -((seq(1.5, 4.0, length.out = 200L) - 2.5)^2) /
          0.05
      )
)
imp_delay <- improve_signal(
  delay_df,
  frequency = 1,
  resample = 1,
  time_min = 0.5,
  time_max = 4.0,
  fourier_components = 0.1
)
expect_true(nrow(imp_delay) > 0L)
expect_true(all(is.finite(imp_delay$intensity)))
expect_true(max(imp_delay$intensity) > 0)
expect_true(all(imp_delay$rtime >= 1.5))

imp_short <- improve_signal(
  data.frame(rtime = c(0.1, 0.12), intensity = c(10, 20)),
  frequency = 1,
  resample = 1,
  time_min = 0.0,
  time_max = 0.15
)
expect_identical(nrow(imp_short), 0L)

# --- preprocess_chromatograms ---
raw_list <- list(
  data.frame(
    rtime = seq(0, 2, length.out = 40L),
    intensity = 5 + 20 * exp(-((seq(0, 2, length.out = 40L) - 1)^2) / 0.05)
  )
)
prep_no_imp <- preprocess_chromatograms(
  detector = "cad",
  list = raw_list,
  name = "sample_A",
  shift = 0.5,
  time_min = 0.2,
  time_max = 1.8,
  baseline_method = "als",
  improve_signal = FALSE
)
expect_identical(
  names(prep_no_imp),
  c(
    "chromatograms_original",
    "chromatograms_original_long",
    "chromatograms_improved",
    "chromatograms_improved_long",
    "chromatograms_baselined",
    "chromatograms_baselined_long"
  )
)
expect_identical(names(prep_no_imp$chromatograms_original), "sample_A")
expect_true(all(
  prep_no_imp$chromatograms_improved$sample_A$rtime >= 0.7 &
    prep_no_imp$chromatograms_improved$sample_A$rtime <= 2.3
))
expect_identical(
  prep_no_imp$chromatograms_original$sample_A$rtime,
  prep_no_imp$chromatograms_original_long$rtime
)
expect_identical(
  prep_no_imp$chromatograms_improved$sample_A$rtime,
  prep_no_imp$chromatograms_improved_long$rtime
)
expect_identical(
  prep_no_imp$chromatograms_baselined$sample_A$rtime,
  prep_no_imp$chromatograms_baselined_long$rtime
)
expect_equal(
  prep_no_imp$chromatograms_original_long$rtime,
  raw_list[[1L]]$rtime + 0.5,
  tolerance = 1e-12
)
expect_equal(min(prep_no_imp$chromatograms_original_long$intensity), 0)
expect_equal(min(prep_no_imp$chromatograms_baselined_long$intensity), 0)
expect_identical(
  prep_no_imp$chromatograms_improved_long$rt_1,
  prep_no_imp$chromatograms_improved_long$rtime
)
expect_identical(
  prep_no_imp$chromatograms_improved_long$rt_2,
  prep_no_imp$chromatograms_improved_long$rtime
)

prep_with_imp <- preprocess_chromatograms(
  detector = "cad",
  list = raw_list,
  name = "sample_A",
  shift = 0.5,
  time_min = 0.2,
  time_max = 1.8,
  fourier_components = 0.2,
  baseline_method = "als",
  improve_signal = TRUE
)
expect_identical(names(prep_with_imp$chromatograms_improved), "sample_A")
expect_true(nrow(prep_with_imp$chromatograms_improved$sample_A) > 0L)

# filter_fft short input (n <= 2L) and high cutoff (components >= 0.5) return x
expect_identical(filter_fft(c(1, 2), components = 0.01), c(1, 2))
expect_identical(filter_fft(c(1, 2, 3, 4), components = 0.5), c(1, 2, 3, 4))
