start <- Sys.time()

pkgload::load_all()

message("This program compares chromatograms.")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

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
