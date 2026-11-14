start <- Sys.time()

out_dir <- file.path(tempdir(), "cascade_peaks")
message("Comparing peaks on bundled example data (writing to ", out_dir, ").")
cascade::process_compare_peaks(
  show_example = TRUE,
  export_dir = out_dir
)

end <- Sys.time()
message("Script finished in ", format(end - start))
