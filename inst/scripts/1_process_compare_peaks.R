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

source(file = "R/change_intensity_name.R")
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
source(file = "R/preprocess_chromatograms.R")
source(file = "R/preprocess_peaks.R")
source(file = "R/process_peaks.R")
source(file = "R/transform_baseline.R")
source(file = "R/transform_ms.R")
source(file = "R/y_as_na.R")

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
  path = file.path(
    paths$inst$extdata$interim$annotations$path,
    params$annotation$tool
  ),
  pattern = params$filename$mzml,
  full.names = TRUE,
  recursive = TRUE
)
FEATURES <- "~/Downloads/test_quant.csv"
EXPORT_DIR <- paths$inst$extdata$interim$peaks
EXPORT_FILE_BPI <-
  paste(params$filename$mzml,
    "featuresInformed",
    "bpi.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_BPI_2 <-
  paste(params$filename$mzml,
    "featuresNotInformed",
    "bpi.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_CAD <-
  paste(params$filename$mzml,
    "featuresInformed",
    "cad.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_CAD_2 <-
  paste(params$filename$mzml,
    "featuresNotInformed",
    "cad.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_PDA <-
  paste(params$filename$mzml,
    "featuresInformed",
    "pda.tsv.gz",
    sep = "_"
  )
EXPORT_FILE_PDA_2 <-
  paste(params$filename$mzml,
    "featuresNotInformed",
    "pda.tsv.gz",
    sep = "_"
  )
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
PEAK_SIMILARITY_PREFILTER <-
  params$chromato$peak$similarity$prefilter
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

log_debug(x = "opening raw files objects and extracting chromatograms")
chromatograms_all <- lapply(files, mzR::openMSfile) |>
  lapply(mzR::chromatograms) |>
  purrr::flatten()

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
data.table::setkey(df_features, rt_1, rt_2)

if (params$signal$detector$bpi == TRUE) {
  detector <- "bpi"
}
if (params$signal$detector$bpi == TRUE) {
  chromatograms_list_bpi <- preprocess_chromatograms(detector = "bpi")
  peaks_prelist_bpi <- preprocess_peaks(detector = "bpi")
  peaks_list_pda <- process_peaks(detector = "bpi")
}
if (params$signal$detector$cad == TRUE) {
  detector <- "cad"
}
if (params$signal$detector$cad == TRUE) {
  chromatograms_list_cad <- preprocess_chromatograms()
  peaks_prelist_cad <- preprocess_peaks()
  peaks_list_cad <- process_peaks()
}
if (params$signal$detector$pda == TRUE) {
  detector <- "pda"
}
if (params$signal$detector$pda == TRUE) {
  #' TODO check if change for PDA doing baselining on normal and not improved
  chromatograms_list_pda <-
    preprocess_chromatograms(detector = "pda")
  peaks_prelist_pda <- preprocess_peaks(detector = "pda")
  peaks_list_pda <- process_peaks(detector = "pda")
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
