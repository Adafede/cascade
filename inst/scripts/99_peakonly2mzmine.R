start <- Sys.time()

library(package = dplyr, quietly = TRUE)
library(package = readr, quietly = TRUE)
source("R/check_export_dir.R")
source("R/get_params.R")
source("R/parse_cli_params.R")
source("R/parse_yaml_params.R")
source("R/parse_yaml_paths.R")

step <- "processing"
paths <- parse_yaml_paths()
params <- ""
params <- get_params(step = step)

log_debug(
  "This program performs",
  "format conversion between peakonlly output and needed MZmine input"
)
log_debug("Authors: \n", "AR")
log_debug("Contributors: \n", "...")

OUTPUT_LIST <-
  list.files(
    path = paths$inst$extdata$interim$peakonly,
    pattern = params$filename$mzml,
    full.names = TRUE
  )

EXPORT <- file.path(
  paths$inst$extdata$interim$mzmine,
  paste(params$filename$mzml, "targeted.csv", sep = "_")
)

list <- lapply(
  OUTPUT_LIST,
  readr::read_delim,
  col_select = c(
    name = 1,
    mz_mean,
    rt_min,
    rt_max,
    intensity = 5
  )
)

table_ready <- list |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    rt_mean = round(x = (rt_min + rt_max) * 30, digits = 2),
    mz = round(x = mz_mean, digits = 2)
  ) |>
  dplyr::filter(intensity >= params$chromato$intensity$ms1$min) |>
  dplyr::distinct(mz,
    rt_mean,
    .keep_all = TRUE
  ) |>
  dplyr::group_by(mz) |>
  dplyr::add_count() |>
  dplyr::ungroup() |>
  dplyr::filter(
    n < params$chromato$occurrence$max &
      rt_mean > params$chromato$time$min &
      rt_mean < params$chromato$time$max
  ) |> ## filter noise
  dplyr::distinct(mz_mean, rt_mean) |>
  dplyr::arrange(mz_mean)

table_ready$name <- row.names(table_ready)

log_debug(x = "checking export directory")
check_export_dir(dirname(EXPORT))
readr::write_delim(
  x = table_ready,
  file = EXPORT,
  delim = ","
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
