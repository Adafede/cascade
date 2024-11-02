start <- Sys.time()

pkgload::load_all()

message("This program compares chromatograms.")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

check_chromatograms_alignment(show_example = TRUE)

end <- Sys.time()

message("Script finished in ", format(end - start))
