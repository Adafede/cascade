start <- Sys.time()

pkgload::load_all()

message("This program generates IDs")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

generate_ids()
# generate_ids(
#   taxa = c("Swertia", "Kopsia", "Ginkgo"),
#   comparison = c("Swertia", "Kopsia"),
#   no_stereo = TRUE,
#   filter_ms_conditions = TRUE,
#   start = "0",
#   end = "9999",
#   limit = "1000000"
#   )

end <- Sys.time()

message("Script finished in ", format(end - start))
