start <- Sys.time()

message("Checking chromatogram alignment on bundled example data.")
cascade::check_chromatograms_alignment(show_example = TRUE)

end <- Sys.time()
message("Script finished in ", format(end - start))
