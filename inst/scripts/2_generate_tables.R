start <- Sys.time()

out_dir <- file.path(tempdir(), "cascade_tables")
if (interactive() || identical(Sys.getenv("CASCADE_RUN_LIVE"), "true")) {
  local({
    old_wd <- setwd(tempdir())
    on.exit(setwd(old_wd), add = TRUE)
    message("Generating summary tables (writing to ", out_dir, ").")
    cascade::generate_tables(show_example = TRUE, export_dir = out_dir)
  })
} else {
  message(
    "Skipping live LOTUS download for generate_tables(); ",
    "run interactively or set CASCADE_RUN_LIVE=true."
  )
}

end <- Sys.time()
message("Script finished in ", format(end - start))
