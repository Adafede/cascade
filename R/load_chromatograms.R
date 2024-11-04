#' Load chromatograms
#'
#' @param file File
#' @param show_example Show example? Default to FALSE
#' @param example_polarity Example polarity
#'
#' @return A list of chromatograms
#'
#' @examples NULL
load_chromatograms <- function(file = NULL,
                               show_example = FALSE,
                               example_polarity = "pos") {
  if (show_example) {
    # chromatograms_positive |>
    #   saveRDS(file = "inst/extdata/chromatograms_positive.rds")
    # chromatograms_negative |>
    #   saveRDS(file = "inst/extdata/chromatograms_negative.rds")
    switch(example_polarity,
      "pos" = readRDS(
        system.file("extdata", "chromatograms_positive.rds", package = "cascade")
      ),
      "neg" = readRDS(
        system.file("extdata", "chromatograms_negative.rds", package = "cascade")
      )
    )
  } else {
    file |>
      mzR::openMSfile() |>
      mzR::chromatograms()
  }
}
