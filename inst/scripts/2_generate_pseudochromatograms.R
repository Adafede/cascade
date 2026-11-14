start <- Sys.time()

message("Generating pseudochromatograms on bundled example data.")
cascade::generate_pseudochromatograms(show_example = TRUE)

end <- Sys.time()
message("Script finished in ", format(end - start))
