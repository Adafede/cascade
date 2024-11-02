start <- Sys.time()

pkgload::load_all()

message("This program compares peaks")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

process_compare_peaks(show_example = TRUE)

end <- Sys.time()

message("Script finished in ", format(end - start))
