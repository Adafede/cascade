start <- Sys.time()

message("Checking peak integration on bundled example data.")
cascade::check_peaks_integration(show_example = TRUE)

end <- Sys.time()
message("Script finished in ", format(end - start))
