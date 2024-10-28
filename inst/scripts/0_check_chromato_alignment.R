start <- Sys.time()

source(file = "R/add_chromato_line.R")
source(file = "R/baseline_chromatogram.R")
source(file = "R/change_intensity_name.R")
source(file = "R/check_chromatograms.R")
source(file = "R/extract_chromatogram.R")
source(file = "R/improve_signals_progress.R")
source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/log_debug.R")
source(file = "R/normalize_chromatograms_list.R")

log_debug("This program compares chromatograms.")
log_debug("Authors: \n", "AR")
log_debug("Contributors: \n", "...")

#' Specific paths
FILE_NEGATIVE <- NULL
FILE_POSITIVE <- "data/source/mzml/210619_AR_06_V_03_2_01.mzML"
TIME_MIN <- 0.7
TIME_MAX <- 35.2
CAD_SHIFT <- 0.05
PDA_SHIFT <- 0.1

if (!is.null(FILE_POSITIVE)) {
  chromatograms_positive <- FILE_POSITIVE |>
    mzR::openMSfile() |>
    mzR::chromatograms()

  chromatogram_bpi_pos <- chromatograms_positive |>
    extract_chromatogram("bpi")
  chromatogram_cad_pos <- chromatograms_positive |>
    extract_chromatogram("cad")
  chromatogram_pda_pos <- chromatograms_positive |>
    extract_chromatogram("pda")

  chromatogram_bpi_pos_improved <- chromatogram_bpi_pos |>
    improve_signal() |>
    normalize_chromatograms_list()
  chromatogram_cad_pos_improved <- chromatogram_cad_pos |>
    improve_signal() |>
    normalize_chromatograms_list()
  chromatogram_pda_pos_improved <- chromatogram_pda_pos |>
    improve_signal() |>
    normalize_chromatograms_list()

  chromatogram_bpi_pos_baselined <- chromatogram_bpi_pos_improved |>
    baseline_chromatogram() |>
    normalize_chromatograms_list()
  chromatogram_cad_pos_baselined <- chromatogram_cad_pos_improved |>
    baseline_chromatogram() |>
    normalize_chromatograms_list()
  chromatogram_pda_pos_baselined <- chromatogram_pda_pos_improved |>
    baseline_chromatogram() |>
    normalize_chromatograms_list()
}

if (!is.null(FILE_NEGATIVE)) {
  chromatograms_negative <- FILE_NEGATIVE |>
    mzR::openMSfile() |>
    mzR::chromatograms()

  chromatogram_bpi_neg <- chromatograms_negative |>
    extract_chromatogram("bpi")
  chromatogram_cad_neg <- chromatograms_negative |>
    extract_chromatogram("cad")
  chromatogram_pda_neg <- chromatograms_negative |>
    extract_chromatogram("pda")

  chromatogram_bpi_neg_improved <- chromatogram_bpi_neg |>
    improve_signal() |>
    normalize_chromatograms_list()
  chromatogram_cad_neg_improved <- chromatogram_cad_neg |>
    improve_signal() |>
    normalize_chromatograms_list()
  chromatogram_pda_neg_improved <- chromatogram_pda_neg |>
    improve_signal() |>
    normalize_chromatograms_list()

  chromatogram_bpi_neg_baselined <- chromatogram_bpi_neg_improved |>
    baseline_chromatogram() |>
    normalize_chromatograms_list()
  chromatogram_cad_neg_baselined <- chromatogram_cad_neg_improved |>
    baseline_chromatogram() |>
    normalize_chromatograms_list()
  chromatogram_pda_neg_baselined <- chromatogram_pda_neg_improved |>
    baseline_chromatogram() |>
    normalize_chromatograms_list()
}

check_chromatograms(
  shift_cad = CAD_SHIFT,
  shift_pda = PDA_SHIFT
)
check_chromatograms(
  shift_cad = CAD_SHIFT,
  shift_pda = PDA_SHIFT,
  type = "baselined"
)
