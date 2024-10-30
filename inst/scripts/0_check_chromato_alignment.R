start <- Sys.time()

pkgload::load_all()

tima::log_debug("This program compares chromatograms.")
tima::log_debug("Authors: \n", "AR")
tima::log_debug("Contributors: \n", "...")

check_chromatograms_alignment(
  show_example = TRUE,
  time_min = 0.5,
  time_max = 35.0,
  type = "improved",
  chromatograms = c(
    "bpi_pos",
    "cad_pos",
    "pda_pos",
    "bpi_neg",
    "cad_neg",
    "pda_neg"
  )
)
