path_mzmine_quant <-
  "~/git/cascade/inst/extdata/interim/mzmine/lists/191109_AR_10043_enriched_UHR_Pos_gnps_quant_full.csv"

path_mzmine_quant_mutated <-
  "~/git/cascade/inst/extdata/interim/mzmine/lists/191109_AR_10043_enriched_UHR_Pos_gnps_quant_full_old.csv"

filename <- path_mzmine_quant |>
  gsub(pattern = "_gnps.*", replacement = "") |>
  gsub(pattern = ".*/", replacement = "") |>
  paste0(".mzML")

fromnewtoold <- function(df) {
  df_new <- df |>
    dplyr::select(
      "row ID" = id,
      "row m/z" = paste0(
        "datafile:",
        filename,
        ":mz"
      ),
      "row retention time" = paste0(
        "datafile:",
        filename,
        ":rt"
      ),
      !!as.name(paste0(
        filename,
        " Peak area"
      )) := paste0(
        "datafile:",
        filename,
        ":area"
      )
    )
}

quant_table <- readr::read_delim(file = path_mzmine_quant)

df_new <- quant_table |>
  fromnewtoold()

readr::write_csv(
  x = df_new,
  file = path_mzmine_quant_mutated
)
