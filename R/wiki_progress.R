#' Wiki progress
#'
#' @param xs XS
#'
#' @return A list of results of Wikidata queries
#'
#' @examples NULL
wiki_progress <- function(xs) {
  xs |>
    purrr::map(
      .f = function(x) {
        WikidataQueryServiceR::query_wikidata(x)
      }
    )
}
