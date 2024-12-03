start <- Sys.time()

pkgload::load_all()

message("This program generates pseudochromatograms")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

generate_pseudochromatograms(show_example = TRUE)
# generate_pseudochromatograms(
#   annotations = "data/interim/annotations/tima/241026_103144_extract/extract_results.tsv",
#   features_informed = "data/interim/peaks/210619_AR_06_V_03_2_01_featuresInformed_cad.tsv",
#   features_not_informed = "data/interim/peaks/210619_AR_06_V_03_2_01_featuresNotInformed_cad.tsv",
#   file = "data/source/mzml/210619_AR_06_V_03_2_01.mzML",
#   headers = c("BasePeak_0", "PDA#1_TotalAbsorbance_0", "UV#1_CAD_1_0"),
#   detector = "cad",
#   show_example = FALSE,
#   min_confidence = 0.4,
#   min_similarity_prefilter = 0.6,
#   min_similarity_filter = 0.8,
#   mode = "pos",
#   organism = "Swertia chirayita",
#   fourier_components = 0.01,
#   frequency = 1,
#   resample = 1,
#   shift = 0.05,
#   time_min = 0.5,
#   time_max = 32.5
#   )

end <- Sys.time()

message("Script finished in ", format(end - start))
