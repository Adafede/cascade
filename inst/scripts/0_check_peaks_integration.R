start <- Sys.time()

pkgload::load_all()

message("This program checks peaks integration.")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

check_peaks_integration(show_example = TRUE)
# check_peaks_integration(
#   file = "data/source/mzml/210619_AR_06_V_03_2_01.mzML",
#   features = "data/interim/mzmine/lists/extract.csv",
#   detector = "cad",
#   chromatogram = "baselined",
#   headers = c("BasePeak_0", "PDA#1_TotalAbsorbance_0", "UV#1_CAD_1_0"),
#   min_area = 0.005,
#   min_intensity = 1E4,
#   shift = 0.05,
#   show_example = FALSE,
#   fourier_components = 0.01,
#   time_min = 0.5,
#   time_max = 32.5,
#   frequency = 1,
#   resample = 1
# )

end <- Sys.time()

message("Script finished in ", format(end - start))
