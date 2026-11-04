library(tinytest)
library(cascade)

peaks_progress <- utils::getFromNamespace("peaks_progress", "cascade")
join_peaks <- utils::getFromNamespace("join_peaks", "cascade")
normalize_chromato <- utils::getFromNamespace("normalize_chromato", "cascade")
prepare_features <- utils::getFromNamespace("prepare_features", "cascade")
prepare_peaks <- utils::getFromNamespace("prepare_peaks", "cascade")
prepare_rt <- utils::getFromNamespace("prepare_rt", "cascade")
prepare_mz <- utils::getFromNamespace("prepare_mz", "cascade")
preprocess_peaks <- utils::getFromNamespace("preprocess_peaks", "cascade")
transform_ms <- utils::getFromNamespace("transform_ms", "cascade")
extract_ms_progress <- utils::getFromNamespace("extract_ms_progress", "cascade")
extract_ms_peak <- utils::getFromNamespace("extract_ms_peak", "cascade")
closest_idx <- utils::getFromNamespace("closest_idx", "cascade")
compare_peaks <- utils::getFromNamespace("compare_peaks", "cascade")

# --- join_peaks boundary tests ---
chrom_dt <- data.table::data.table(
  id = rep("s1", 6L),
  rt_1 = c(1, 2, 3, 5, 6, 7),
  rt_2 = c(1, 2, 3, 5, 6, 7),
  intensity = c(10, 20, 10, 30, 60, 30)
)
peaks_dt <- data.table::data.table(
  id = c("s1", "s1"),
  peak_id = c(1L, 2L),
  peak_max = c(20, 60),
  rt_apex = c(2, 6),
  rt_min = c(1, 5),
  rt_max = c(3, 7)
)
# Peak 1 integral = 40, Peak 2 integral = 120, total = 160 -> fractions 0.25 and 0.75
jp_exact <- join_peaks(chrom_dt, peaks_dt, min_area = 0.25)
expect_identical(nrow(jp_exact), 2L)
expect_identical(jp_exact$integral, c(40, 120))

jp_above <- join_peaks(chrom_dt, peaks_dt, min_area = 0.2501)
expect_identical(nrow(jp_above), 1L)
expect_identical(jp_above$peak_id, 2L)
expect_identical(jp_above$integral, 120)

# --- normalize_chromato boundary tests ---
df_xy <- data.frame(
  rtime = c(1, 2, 3, 4, 5),
  intensity = c(10, 20, 60, 110, 10)
)
# Slice [2, 4] has intensities c(20, 60, 110) -> normalized c(0, 40/90, 1)
pk_bounds <- data.frame(rt_min = 2, rt_max = 4)
nc_res <- normalize_chromato(pk_bounds, df_xy, intensity_threshold = 40 / 90)
expect_identical(nc_res$rtime, c(3, 4))
expect_equal(nc_res$intensity, c(40 / 90, 1.0), tolerance = 1e-12)

nc_above <- normalize_chromato(pk_bounds, df_xy, intensity_threshold = 0.5)
expect_identical(nc_above$rtime, 4)
expect_equal(nc_above$intensity, 1.0, tolerance = 1e-12)

# Empty slice
nc_empty <- normalize_chromato(data.frame(rt_min = 10, rt_max = 12), df_xy)
expect_identical(nrow(nc_empty), 0L)
expect_identical(names(nc_empty), c("rtime", "intensity"))

# --- prepare_features ---
sample_nm <- "sample.mzML"
raw_feat <- data.frame(
  id = c("101", "102", "103"),
  rt = c("1.5", "2.5", "3.5"),
  mz = c("300.1", "400.2", "500.3"),
  area = c("1000", "2000", "3000"),
  status = c("DETECTED", "ESTIMATED", "DETECTED"),
  rt_min = c("1.4", "2.4", "3.4"),
  rt_max = c("1.6", "2.6", "3.6"),
  mz_min = c("300.09", "400.19", "500.29"),
  mz_max = c("300.11", "400.21", "500.31"),
  int_min = c("10", "20", "30"),
  int_max = c("500", "9000", "1000"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(raw_feat)[5:11] <- c(
  paste0("datafile:", sample_nm, ":feature_state"),
  paste0("datafile:", sample_nm, ":rt_range:min"),
  paste0("datafile:", sample_nm, ":rt_range:max"),
  paste0("datafile:", sample_nm, ":mz_range:min"),
  paste0("datafile:", sample_nm, ":mz_range:max"),
  paste0("datafile:", sample_nm, ":intensity_range:min"),
  paste0("datafile:", sample_nm, ":intensity_range:max")
)
pf_500 <- prepare_features(raw_feat, min_intensity = 500, name = sample_nm)
expect_identical(pf_500$feature_id, c(101, 103))
expect_identical(data.table::key(pf_500), c("rt_1", "rt_2"))

pf_501 <- prepare_features(raw_feat, min_intensity = 501, name = sample_nm)
expect_identical(pf_501$feature_id, 103)

# --- prepare_peaks, prepare_rt, prepare_mz ---
pp_res <- prepare_peaks(data.frame(rtime = c(2, 4, 10), intensity = c(5, 9, 7)))
expect_equal(pp_res$rtime, c(0, 0.25, 1), tolerance = 1e-12)
expect_identical(pp_res$intensity, c(5, 9, 7))

pp_empty <- prepare_peaks(data.frame(
  rtime = numeric(0L),
  intensity = numeric(0L)
))
expect_identical(pp_empty$rtime, numeric(0L))
expect_identical(pp_empty$intensity, numeric(0L))

prt <- prepare_rt(
  data.frame(rt_min = c(1.5, 9), rt_max = c(2.5, 9)),
  shift = 0.5
)
expect_identical(prt, cbind(rtmin = 120, rtmax = 180))

pmz <- prepare_mz(data.frame(
  mz_min = c(100.1, 200.2),
  mz_max = c(100.3, 200.4)
))
expect_identical(
  pmz,
  list(
    cbind(mzmin = 100.1, mzmax = 100.3),
    cbind(mzmin = 200.2, mzmax = 200.4)
  )
)

# --- preprocess_peaks with 0 overlapping features vs overlapping features ---
t_seq <- seq(0, 10, by = 0.1)
y_seq <- 100 * exp(-((t_seq - 5)^2) / 0.2)
df_xy_synth <- data.frame(rtime = t_seq, intensity = y_seq)
df_long_synth <- data.table::data.table(
  id = "1",
  rtime = t_seq,
  intensity = y_seq,
  rt_1 = t_seq,
  rt_2 = t_seq
)
# Non-overlapping feature at rt [20.0, 25.0] (outside chromatogram [0, 10])
feat_no_ov <- data.table::data.table(
  feature_id = 1,
  rt = 22.0,
  mz = 250.1,
  area = 500,
  rt_1 = 20.0,
  rt_2 = 25.0,
  mz_min = 250.09,
  mz_max = 250.11,
  intensity_min = 10,
  intensity_max = 1000
)
data.table::setkey(feat_no_ov, rt_1, rt_2)
pp_no_ov <- preprocess_peaks(
  detector = "cad",
  df_features = feat_no_ov,
  df_long = df_long_synth,
  df_xy = df_xy_synth,
  name = "sample.mzML",
  shift = 0,
  min_area = 0.01,
  fit = "gaussian"
)
expect_identical(length(pp_no_ov$list_df_features_with_peaks_long), 0L)
expect_identical(nrow(pp_no_ov$df_features_without_peaks), 1L)

# --- transform_ms & extract_ms_peak ---
raw_xic <- list(
  list(
    intensity = c(NA_real_, 10, 30, 50, 110),
    rtime = c(1, 5, 2, 4, 3)
  ),
  list(
    intensity = numeric(0L),
    rtime = numeric(0L)
  )
)
tms <- transform_ms(raw_xic, min_intensity = 0.2)
# Valid intensities: 10, 30, 50, 110 -> norm: 0, 0.2, 0.4, 1.0
# Keep >= 0.2 -> intensities 0.2, 0.4, 1.0 at rtime 2, 4, 3
# Norm rtime of (2, 4, 3) -> (0, 1, 0.5); sorted by rtime -> rtimes (0, 0.5, 1), intensities (0.2, 1.0, 0.4)
expect_equal(tms[[1L]]$rtime, c(0, 0.5, 1), tolerance = 1e-12)
expect_equal(tms[[1L]]$intensity, c(0.2, 1.0, 0.4), tolerance = 1e-12)
expect_identical(nrow(tms[[2L]]), 0L)

emp <- extract_ms_peak(tms)
expect_equal(emp[[1L]]$rtime, c(0, 0.5, 1), tolerance = 1e-12)
expect_equal(emp[[1L]]$intensity, c(0.2, 1.0, 0.4), tolerance = 1e-12)

# --- extract_ms_progress branches (empty RT window, empty scan, 1 hit, >1 hits) ---
ms_synth <- list(
  rt = c(10, 20, 30, 40),
  mz = list(
    numeric(0L),
    c(100.0, 200.0, 200.05),
    c(100.0, 200.02),
    c(100.0, 200.01, 200.04)
  ),
  intensity = list(
    numeric(0L),
    c(10, 100, 50),
    c(30, 200),
    c(50, 300, 100)
  )
)
rts_list <- list(
  cbind(rtmin = 5, rtmax = 8), # Peak 1: no scans in window (r_rt <= l_rt)
  cbind(rtmin = 10, rtmax = 40) # Peak 2: scans 1..4 (scan 1 is empty, scans 2..4 have 1 or 2 hits)
)
mzs_list <- list(
  list(cbind(mzmin = 99.9, mzmax = 100.1)),
  list(
    cbind(mzmin = 99.9, mzmax = 100.1), # Feature 1: single hit per scan (10, 30, 50)
    cbind(mzmin = 199.9, mzmax = 200.1) # Feature 2: multi-hit summed per scan (150, 200, 400)
  )
)
ext_res <- extract_ms_progress(
  xs = 1:2,
  ms_data = ms_synth,
  rts = rts_list,
  mzs = mzs_list
)
expect_identical(nrow(ext_res[[1L]][[1L]]), 0L)
# Feature 1 in Peak 2: raw intensities (NA, 10, 30, 50) at rt (10, 20, 30, 40) -> norm int (0, 0.5, 1) -> keep >= 0.1 gives (0.5, 1) at rt (0, 1)
expect_equal(ext_res[[2L]][[1L]]$intensity, c(0.5, 1.0), tolerance = 1e-12)
expect_equal(ext_res[[2L]][[1L]]$rtime, c(0, 1), tolerance = 1e-12)
# Feature 2 in Peak 2: raw summed intensities (NA, 150, 200, 400) -> norm int (0, 50/250=0.2, 1.0) -> keep >= 0.1 gives (0.2, 1.0)
expect_equal(ext_res[[2L]][[2L]]$intensity, c(0.2, 1.0), tolerance = 1e-12)
expect_equal(ext_res[[2L]][[2L]]$rtime, c(0, 1), tolerance = 1e-12)

# --- closest_idx all branches ---
expect_identical(closest_idx(numeric(0L), c(1, 2), 0.5), integer(0L))
expect_identical(
  closest_idx(c(1, 2), numeric(0L), 0.5),
  c(NA_integer_, NA_integer_)
)
# Queries remaining after table exhausted (it > nt)
expect_identical(
  closest_idx(c(1.0, 5.0, 9.0), c(1.0), tolerance = 0.2),
  c(1L, NA_integer_, NA_integer_)
)
# Duplicate match to same table element: second x is closer so first is overwritten with NA
expect_identical(
  closest_idx(c(0.9, 0.99), c(1.0, 5.0), tolerance = 0.2),
  c(NA_integer_, 1L)
)
# Duplicate match to same table element: first x is closer than second x
expect_identical(
  closest_idx(c(1.01, 1.15), c(1.0, 5.0), tolerance = 0.2),
  c(1L, NA_integer_)
)
# Exact tolerance boundary (including IEEE 754 positive rounding 0.3 - 0.1 > 0.2)
expect_identical(
  closest_idx(c(0.3, 0.35), c(0.1), tolerance = 0.2),
  c(1L, NA_integer_)
)

# --- compare_peaks all branches ---
# 1. Empty ms_pks
expect_identical(
  compare_peaks(1L, list(list()), list()),
  list(0)
)
# 2. Identical retention time vectors (direct correlation without closest_idx)
p_pre_same <- list(
  list_chromato_peaks = list(
    list(rtime = c(0, 0.5, 1), intensity = c(1, 2, 3))
  )
)
ms_same <- list(
  list(
    list(rtime = c(0, 0.5, 1), intensity = c(2, 4, 6)),
    list(rtime = c(0, 0.5, 1), intensity = c(6, 4, 2))
  )
)
cp_same <- compare_peaks(1L, ms_same, p_pre_same)
expect_equal(cp_same[[1L]], 1.0, tolerance = 1e-12)
expect_equal(cp_same[[2L]], -1.0, tolerance = 1e-12)

# 3. Different retention time vectors aligned via closest_idx
p_pre_diff <- list(
  list_chromato_peaks = list(
    list(rtime = c(0.0, 0.25, 0.5, 0.75, 1.0), intensity = c(1, 2, 4, 2, 1))
  )
)
ms_diff <- list(
  list(
    list(rtime = c(0.01, 0.51, 0.99), intensity = c(10, 40, 10))
  )
)
cp_diff <- compare_peaks(1L, ms_diff, p_pre_diff)
expect_equal(cp_diff[[1L]], 1.0, tolerance = 1e-12)
