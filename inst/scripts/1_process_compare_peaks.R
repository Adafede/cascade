start <- Sys.time()

pkgload::load_all()

message("This program compares peaks")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

process_compare_peaks(show_example = TRUE)
# process_compare_peaks(
#   file = "data/source/mzml/210619_AR_06_V_03_2_01.mzML",
#   features = "data/interim/mzmine/lists/extract.csv",
#   type = "baselined",
#   detector = "cad",
#   headers = c("BasePeak_0", "PDA#1_TotalAbsorbance_0", "UV#1_CAD_1_0"),
#   export_dir = "data/interim/peaks",
#   show_example = FALSE,
#   fourier_components = 0.01,
#   frequency = 1,
#   min_area = 0.005,
#   min_intensity = 1E4,
#   resample = 1,
#   shift = 0.05,
#   time_min = 0.5,
#   time_max = 32.5
#   )

end <- Sys.time()

message("Script finished in ", format(end - start))
