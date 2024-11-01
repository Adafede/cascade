start <- Sys.time()

pkgload::load_all()

message("This program checks peaks integration.")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

check_peaks_integration(show_example = TRUE)

end <- Sys.time()

message("Script finished in ", format(end - start))
