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
    dataset_name <- switch(
      example_polarity,
      "pos" = "cascade_chromatograms_positive",
      "neg" = "cascade_chromatograms_negative"
    )
    utils::data(list = dataset_name, package = "cascade", envir = environment())
    chromatograms_df <- get(dataset_name)
    chromatograms <- lapply(
      names(headers),
      function(chrom_name) {
        chromatograms_df |>
          tidytable::filter(chromatogram == chrom_name) |>
          tidytable::select(rtime, intensity) |>
          data.frame()
      }
    )
    names(chromatograms) <- names(headers)
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
