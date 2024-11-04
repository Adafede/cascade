start <- Sys.time()

pkgload::load_all()

message("This program plots IDs")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

generate_ids()

end <- Sys.time()

message("Script finished in ", format(end - start))
