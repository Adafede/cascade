start <- Sys.time()

if (interactive() || identical(Sys.getenv("CASCADE_RUN_LIVE"), "true")) {
  local({
    old_wd <- setwd(tempdir())
    on.exit(setwd(old_wd), add = TRUE)
    message(
      "Generating literature chemotaxonomic hierarchies via Wikidata/LOTUS."
    )
    cascade::generate_ids()
  })
} else {
  message(
    "Skipping live Wikidata/Zenodo query; ",
    "run interactively or set CASCADE_RUN_LIVE=true."
  )
}

end <- Sys.time()
message("Script finished in ", format(end - start))
