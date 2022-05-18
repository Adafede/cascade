require(dplyr)

#' Title
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
no_other <- function(dataframe) {
  dataframe_no_other <- dataframe |>
    dplyr::filter(!grepl(
      pattern = "not",
      x = best_candidate_1
    ) &
      !is.na(best_candidate_1))

  return(dataframe_no_other)
}
