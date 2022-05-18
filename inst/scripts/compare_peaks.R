start <- Sys.time()

library(package = baseline, quietly = TRUE)
library(package = data.table, quietly = TRUE)
library(package = docopt, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = future, quietly = TRUE)
library(package = future.apply, quietly = TRUE)
library(package = MSnbase, quietly = TRUE)
library(package = nucleR, quietly = TRUE)
library(package = parallel, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = pracma, quietly = TRUE)
library(package = progressr, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = xcms, quietly = TRUE)

source(file = "R/check_export_dir.R")
source(file = "R/compare_peaks.R")
source(file = "R/extract_ms.R")
source(file = "R/extract_ms_peak.R")
source(file = "R/filter_ms.R")
source(file = "R/get_gnps.R")
source(file = "R/get_params.R")
source(file = "R/improve_signal.R")
source(file = "R/improve_signals_progress.R")
source(file = "R/join_peaks.R")
source(file = "R/log_debug.R")
source(file = "R/make_confident.R")
source(file = "R/normalize_chromato.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/peaks_progress.R")
source(file = "R/plot_chromatogram.R")
source(file = "R/prepare_mz.R")
source(file = "R/prepare_peaks.R")
source(file = "R/prepare_rt.R")
source(file = "R/transform_baseline.R")
source(file = "R/transform_ms.R")
source(file = "R/y_as_na.R")

future::plan(strategy = future::multisession)
progressr::handlers(global = TRUE)
progressr::handlers("progress")

step <- "processing"
paths <- parse_yaml_paths()
params <- ""
params <- get_params(step = step)

log_debug(
  "This program performs",
  "Quantitative and Qualitative Contextualization",
  "of in depth annotated extracts"
)
log_debug("Authors: \n", "AR")
log_debug("Contributors: \n", "...")

#' Paths
ANNOTATIONS <- list.files(
  path = file.path(paths$inst$extdata$interim$annotations$path, params$annotation$tool),
  pattern = params$filename$mzml,
  full.names = TRUE,
  recursive = TRUE
)
FEATURES <- "~/Downloads/test_quant.csv"
EXPORT_DIR <- paths$inst$extdata$interim$peaks
EXPORT_FILE <- paste(params$filename$mzml, "featuresInformed.tsv.gz", sep = "_")
EXPORT_FILE_2 <- paste(params$filename$mzml, "featuresNotInformed.tsv.gz", sep = "_")
# GNPS_JOB <- "97d7c50031a84b9ba2747e883a5907cd"
# TOYSET <- "~/data/20210701_10043/fractions"
# TOYSET <- "~/../../Volumes/LaCie/data/20210701_10043/fractions"
# TOYSET <- "~/data/20210701_10043/test"
TOYSET <- paths$inst$extdata$source$mzml$path

#' Generic parameters
WORKERS <- params$workers

#' Parameters for LC alignment
TIME_MIN <- params$chromato$time$min
TIME_MAX <- params$chromato$time$max
CAD_SHIFT <- params$chromato$shift$cad
PDA_SHIFT <- params$chromato$shift$pda
ESTIMATED_SOLUBLITIY_LIMIT <- params$misc$solubility$limit

#' Parameters for signal improvement
FOURRIER_COMPONENTS <- params$signal$fourrier$components
FREQUENCY <- params$signal$frequency
RESAMPLE <- params$signal$resample

#' Parameters adapted from Excel sheet from paper shortDOI: 10/ghmvhz
sigma <- params$signal$sigma
k2 <- sigma / params$signal$k2 # 30
k4 <- sigma / params$signal$k4 # 200
smoothing_width <- params$signal$smoothing
baseline_adjust <- params$signal$baseline

#' Parameters related to MS/CAD
INTENSITY_MS_MIN <- params$chromato$intensity$ms1$min
PEAK_SIMILARITY <- params$chromato$peak$similarity$filter
PEAK_SIMILARITY_PREFILTER <- params$chromato$peak$similarity$prefilter
RT_TOL <- params$chromato$peak$tolerance$rt
PPM <- params$chromato$peak$tolerance$ppm
AREA_MIN <- params$chromato$peak$area$min

#' Parameters for annotation
CONFIDENCE_SCORE_MIN <- params$annotation$confidence$min

log_debug(x = "listing files")
files <- list.files(
  path = TOYSET,
  pattern = ".mzML",
  full.names = TRUE,
  recursive = TRUE
)

# files <- files[grepl(pattern = "M_17|M_28|M_36|M_40|M_47|M_57|M_67", x = files)]

names <- list.files(
  path = TOYSET,
  pattern = ".mzML",
  recursive = TRUE
) |>
  gsub(pattern = "[0-9]{6}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML",
    replacement = "",
    fixed = TRUE
  )

# names <- names[grepl(pattern = "M_17|M_28|M_36|M_40|M_47|M_57|M_67", x = names)]

log_debug(x = "loading feature table")
# feature_table <- read_features(id = GNPS_JOB)
feature_table <- readr::read_delim(file = FEATURES)

log_debug(x = "loading raw files (can take long if loading multiple files)")
dda_data <- MSnbase::readMSData(
  files = files,
  mode = "onDisk",
  msLevel. = 1
)

log_debug(x = "opening raw files objects")
objects <- lapply(files, mzR::openMSfile)

log_debug(x = "extracting chromatograms")
chromatograms <- lapply(objects, mzR::chromatograms)

log_debug(x = "selecting specific chromatograms ...")
chromatograms_all <- chromatograms |>
  purrr::flatten()
log_debug(x = "... CAD")
chromatograms_cad <- chromatograms_all[c(FALSE, FALSE, TRUE)]
log_debug(x = "... PDA")
chromatograms_pda <- chromatograms_all[c(FALSE, TRUE, FALSE)]

log_debug(x = "preparing chromatograms ...")
log_debug(x = "... CAD")
chromatograms_cad_ready <- lapply(chromatograms_cad, function(x) {
  x |>
    dplyr::select(time,
      intensity = UV.1_CAD_1_0
    )
})
log_debug(x = "... PDA")
chromatograms_pda_ready <- lapply(chromatograms_pda, function(x) {
  x |>
    dplyr::select(time,
      intensity = PDA.1_TotalAbsorbance_0
    )
})

log_debug(x = "improving chromatograms ...")
log_debug(x = "... CAD")
chromatograms_cad_improved <-
  improve_signals_progress(chromatograms_cad_ready)
log_debug(x = "... PDA")
chromatograms_pda_improved <-
  improve_signals_progress(chromatograms_pda_ready)
chromatograms_pda_improved <- chromatograms_pda_ready

names(chromatograms_cad_improved) <- names
names(chromatograms_pda_improved) <- names

cads_improved <-
  dplyr::bind_rows(chromatograms_cad_improved, .id = "id") |>
  dplyr::mutate(time = time + CAD_SHIFT)
pdas_improved <-
  dplyr::bind_rows(chromatograms_pda_improved, .id = "id") |>
  dplyr::mutate(time = time + PDA_SHIFT)

log_debug(x = "plotting improved chromatograms ...")
log_debug(x = "... CAD")
cad_plot <- plot_chromatogram(df = cads_improved, text = "CAD")
log_debug(x = "... PDA")
pda_plot <- plot_chromatogram(df = pdas_improved, text = "PDA")

log_debug(x = "baselining chromatograms ...")
log_debug(x = "... CAD")
chromatograms_cad_baselined <-
  transform_baseline(x = chromatograms_cad_improved)
log_debug(x = "... PDA")
chromatograms_pda_baselined <-
  transform_baseline(x = chromatograms_pda_improved)

cads_baselined <-
  dplyr::bind_rows(chromatograms_cad_baselined, .id = "id") |>
  dplyr::mutate(intensity = intensity / max(intensity))
pdas_baselined <-
  dplyr::bind_rows(chromatograms_pda_baselined, .id = "id") |>
  dplyr::mutate(intensity = intensity / max(intensity))

log_debug(x = "plotting baselined chromatograms ...")
log_debug(x = "... CAD")
new_new_cad <- plot_chromatogram(df = cads_baselined, text = "CAD")
log_debug(x = "... PDA")
new_new_pda <- plot_chromatogram(df = pdas_baselined, text = "PDA")

#' data.table call outside of future because buggy else
peaks_cad <- peaks_progress(chromatograms_cad_baselined)

peaks_pda <- peaks_progress(chromatograms_pda_baselined)

names(peaks_cad) <- names
names(peaks_pda) <- names

#' data.table call outside of future because buggy else
peaks_cad_all <- dplyr::bind_rows(peaks_cad, .id = "id") |>
  data.table::data.table()
peaks_pda_all <- dplyr::bind_rows(peaks_pda, .id = "id") |>
  data.table::data.table()

cads_baselined <- cads_baselined |>
  dplyr::mutate(rt_1 = time, rt_2 = time) |>
  data.table::data.table()
pdas_baselined <- pdas_baselined |>
  dplyr::mutate(rt_1 = time, rt_2 = time) |>
  data.table::data.table()

log_debug(x = "joining peaks ...")
log_debug(x = "... CAD")
df_cad <-
  join_peaks(chromatograms = cads_baselined, peaks = peaks_cad_all)
log_debug(x = "... PDA")
df_pda <-
  join_peaks(chromatograms = pdas_baselined, peaks = peaks_pda_all)

log_debug(x = "preparing feature list ...")
colnames(feature_table) <-
  gsub(
    pattern = ".Peak.area",
    replacement = "",
    x = colnames(feature_table)
  )

log_debug(x = "... removing \"row m/z\" and from \"row retention time\" columns")
feature_table <- feature_table |>
  dplyr::select(
    rt = "row retention time",
    mz = "row m/z",
    everything(),
    -"correlation group ID",
    -"annotation network number",
    -"best ion",
    -"auto MS2 verify",
    -"identified by n=",
    -"partners",
    -"neutral M mass"
  ) |>
  dplyr::select(-(ncol(feature_table) - 7)) |>
  tibble::column_to_rownames(var = "row ID")

log_debug(x = "... keeping features above desired intensity only")
df_features <- feature_table |>
  tibble::rownames_to_column() |>
  dplyr::group_by(rowname, rt, mz) |>
  tidyr::gather(column, value, -rowname, -rt, -mz) |>
  dplyr::ungroup() |>
  dplyr::mutate(column = gsub(
    pattern = "^X",
    replacement = "",
    x = column
  )) |>
  dplyr::arrange(rowname, dplyr::desc(value)) |>
  dplyr::filter(value >= INTENSITY_MS_MIN) |>
  dplyr::mutate(column = gsub(
    pattern = ".Peak.area",
    replacement = "",
    x = column
  )) |>
  dplyr::select(
    feature_id = rowname,
    rt,
    mz,
    sample = column,
    intensity = value
  ) |>
  dplyr::mutate(
    rt_1 = as.numeric(rt),
    rt_2 = as.numeric(rt),
    mz_min = (1 - (1E-6 * PPM)) * as.numeric(mz),
    mz_max = (1 + (1E-6 * PPM)) * as.numeric(mz)
  ) |>
  data.table::data.table()

log_debug(x = "setting joining keys")
data.table::setkey(df_cad, rt_min, rt_max)
data.table::setkey(df_features, rt_1, rt_2)

log_debug(x = "joining within given rt tolerance")
df_new_pre <- data.table::foverlaps(df_features, df_cad)

log_debug(x = "selecting features within peaks")
df_new_with <- df_new_pre |>
  dplyr::select(-rt_1, -rt_2) |>
  dplyr::filter(!is.na(peak_id))

# df_new_with <- df_new_with |> #' TODO DONT FORGET
#   dplyr::distinct(id, peak_id, feature_id, .keep_all = TRUE) |> #' TODO DONT FORGET
#   sample_n(500) #' TODO DONT FORGET

log_debug(x = "selecting features outside peaks")
df_new_without <- df_new_pre |>
  dplyr::filter(is.na(peak_id)) |>
  dplyr::filter(sample %in% df_new_with$sample)

# df_new_without <- df_new_without |> #' TODO DONT FORGET
#   dplyr::distinct(id, peak_id, feature_id, .keep_all = TRUE) #' TODO DONT FORGET

log_debug(x = "splitting by file")
list_df_peaks <- df_new_with |>
  dplyr::group_split(id)

names(list_df_peaks) <- unique(df_new_with$id)

log_debug(x = "splitting by peak")
list_df_peaks_split <- parallel::mclapply(
  X = list_df_peaks,
  FUN = dplyr::group_split,
  peak_id
)

list_df_peaks_split_flat <- list_df_peaks_split |>
  purrr::flatten()

log_debug(x = "retrieving ms files") #' TODO warning CAD shift only
list_dda_split <- parallel::mclapply(
  X = list_df_peaks_split_flat,
  FUN = filter_ms
)

log_debug(x = "normalizing chromato")
list_chromato_split <-
  parallel::mclapply(
    X = list_df_peaks_split_flat,
    FUN = normalize_chromato
  )

log_debug(x = "preparing peaks chromato")
list_chromato_peaks <- parallel::mclapply(
  X = list_chromato_split,
  FUN = prepare_peaks
)

log_debug(x = "preparing rt")
list_rtr <- parallel::mclapply(
  X = list_df_peaks_split_flat,
  FUN = prepare_rt
)

log_debug(x = "preparing mz")
list_mzr <- parallel::mclapply(
  X = list_df_peaks_split_flat,
  FUN = prepare_mz
)

log_debug(x = "extracting ms chromatograms (longest step)")
log_debug(x = "count approx 1 minute per 100 features (increasing with features number)")
log_debug(x = "varies a lot depending on features distribution")
list_ms_chr <-
  parallel::mclapply(
    X = seq_along(list_df_peaks_split_flat),
    FUN = extract_ms
  )

log_debug(x = "transforming ms chromatograms")
list_ms_transformed <- parallel::mclapply(
  X = list_ms_chr,
  FUN = transform_ms
)

log_debug(x = "extracting ms peaks")
list_ms_peaks <- parallel::mclapply(
  X = list_ms_transformed,
  FUN = extract_ms_peak
)

log_debug(x = "comparing peaks")
list_comparison_score <-
  parallel::mclapply(
    X = seq_along(list_ms_peaks),
    FUN = compare_peaks
  )

df_peaks_samples_full <- list_df_peaks_split_flat |>
  dplyr::bind_rows()

comparison_scores <- list_comparison_score |>
  purrr::flatten()

df_peaks_samples_full$comparison_score <-
  as.numeric(comparison_scores)

log_debug(x = "final aesthetics")
df_final <- df_peaks_samples_full |>
  dplyr::select(
    sample = id,
    peak_id,
    peak_rt_min = rt_min,
    peak_rt_apex = rt_apex,
    peak_rt_max = rt_max,
    peak_area = integral,
    feature_id,
    feature_rt = rt,
    feature_mz = mz,
    feature_area = intensity,
    comparison_score
  ) |>
  dplyr::distinct()

df_final_2 <- df_new_without |>
  dplyr::mutate(comparison_score = NA) |>
  dplyr::select(
    sample = id,
    peak_id,
    peak_rt_min = rt_min,
    peak_rt_apex = rt_apex,
    peak_rt_max = rt_max,
    peak_area = integral,
    feature_id,
    feature_rt = rt,
    feature_mz = mz,
    feature_area = intensity,
    comparison_score
  ) |>
  dplyr::distinct()

log_debug(x = "checking export directory")
check_export_dir(EXPORT_DIR)

log_debug(x = "exporting to ...")
log_debug(x = EXPORT_DIR)
readr::write_tsv(
  x = df_final,
  file = file.path(EXPORT_DIR, EXPORT_FILE)
)
readr::write_tsv(
  x = df_final_2,
  file = file.path(EXPORT_DIR, EXPORT_FILE_2)
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
