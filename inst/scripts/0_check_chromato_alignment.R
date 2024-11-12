start <- Sys.time()

pkgload::load_all()

message("This program compares chromatograms.")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

check_chromatograms_alignment(show_example = TRUE)
# check_chromatograms_alignment(
#   file_negative = NULL,
#   file_positive = "data/source/mzml/210619_AR_06_V_03_2_01.mzML",
#   time_min = 0.5,
#   time_max = 32.5,
#   cad_shift = 0.05,
#   pda_shift = 0.1,
#   fourier_components = 0.01,
#   frequency = 1,
#   resample = 1,
#   chromatograms = c("bpi_pos", "cad_pos", "pda_pos"),
#   type = "baselined",
#   normalize_intensity = TRUE,
#   normalize_time = FALSE,
#   show_example = FALSE
# )

end <- Sys.time()

message("Script finished in ", format(end - start))
