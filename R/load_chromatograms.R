#' Load chromatograms
#'
#' @param file File
#' @param headers Headers
#' @param show_example Show example? Default to FALSE
#' @param example_polarity Example polarity
#'
#' @return A list of chromatograms
#'
#' @examples NULL
load_chromatograms <- function(
  file = NULL,
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  show_example = FALSE,
  example_polarity = "pos"
) {
  if (show_example) {
    # chromatograms_positive |>
    #   saveRDS(file = "inst/extdata/chromatograms_positive.rds")
    # chromatograms_negative |>
    #   saveRDS(file = "inst/extdata/chromatograms_negative.rds")
    chromatograms <- switch(
      example_polarity,
      "pos" = readRDS(
        system.file(
          "extdata",
          "chromatograms_positive.rds",
          package = "cascade"
        )
      ),
      "neg" = readRDS(
        system.file(
          "extdata",
          "chromatograms_negative.rds",
          package = "cascade"
        )
      )
    )
    names(chromatograms) <- headers |>
      names()
    chromatograms <- chromatograms |>
      ## TODO dirty fix
      purrr::map(
        .f = function(df) {
          df |> tidytable::rename(rtime = 1, intensity = 2)
        }
      )
    return(chromatograms)
  } else {
    file_pointer <- file |>
      mzR::openMSfile()
    file_headers <- file_pointer |>
      mzR::chromatogramHeader()
    indices <- file_headers$chromatogramIndex[
      file_headers$chromatogramId %in% headers
    ]
    chromatograms <- file_pointer |>
      mzR::chromatograms(chrom = indices)
    names(chromatograms) <- headers |>
      names()
    return(chromatograms)
  }
}
