#' No other
#'
#' @param dataframe Dataframe
#'
#' @return A dataframe with no other
#'
#' @examples NULL
no_other <- function(dataframe) {
  dataframe |>
    tidytable::filter(!grepl(
      pattern = "not",
      x = best_candidate_1
    ) &
      !is.na(best_candidate_1))
}
