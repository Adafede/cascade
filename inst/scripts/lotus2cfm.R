source("paths.R")

export_path <- "~/Downloads/tmp/lotus/smiles4cfm.txt"

message("exporting unique LOTUS 2D structures for yaccl classification")
lotus <- readr::read_delim(
  file = file.path(
    pathDataProcessed,
    pathLastFrozen
  ),
  col_types = cols(.default = "c"),
  locale = locales,
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
