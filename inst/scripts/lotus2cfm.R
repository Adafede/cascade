start <- Sys.time()

#' Packages
packages_cran <-
  c(
    "dplyr",
    "readr",
    "stringr",
    "yaml"
  )
packages_bioconductor <- NULL
packages_github <- NULL

source(file = "R/check_and_load_packages.R")
source(file = "R/load_lotus.R")
source(file = "R/parse_yaml_paths.R")

paths <- parse_yaml_paths()

export_path <- "~/Downloads/tmp/lotus/smiles4cfm.txt"

load_lotus()

message("exporting unique LOTUS 2D structures for yaccl classification")
lotus <- readr::read_delim(
  file = paths$inst$extdata$source$libraries$lotus,
  col_select = c("structure_inchikey", "structure_smiles_2D")
) |>
  dplyr::distinct() |>
  dplyr::mutate(structure_inchikey_2D = gsub(
    pattern = "-.*",
    replacement = "",
    x = structure_inchikey
  )) |>
  dplyr::select(structure_inchikey_2D, structure_smiles_2D) |>
  dplyr::distinct()

readr::write_delim(
  x = lotus,
  file = export_path,
  col_names = FALSE
)
