start <- Sys.time()

pkgload::load_all()

message("This program plots pseudochromatograms")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

process_plot_pseudochromatograms(show_example = TRUE)

end <- Sys.time()

message("Script finished in ", format(end - start))
