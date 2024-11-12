start <- Sys.time()

pkgload::load_all()

message("This program generates tables")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

generate_tables(show_example = TRUE)
# generate_tables(
#   annotations = "data/interim/annotations/tima/241026_103144_extract/extract_results.tsv",
#   file_negative = NULL,
#   file_positive = "data/interim/peaks/210619_AR_06_V_03_2_01_featuresInformed_cad.tsv",
#   min_confidence = 0.4,
#   show_example = FALSE,
#   export_csv = TRUE,
#   export_html = TRUE,
#   export_dir = "data/processed",
#   export_name = "cascade_table"
# )

end <- Sys.time()

message("Script finished in ", format(end - start))
