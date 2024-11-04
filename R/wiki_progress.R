#' Wiki progress
#'
#' @param xs XS
#'
#' @return A list of results of Wikidata queries
#'
#' @examples NULL
wiki_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  lapply(
    X = xs,
    FUN = function(x) {
      p()
      WikidataQueryServiceR::query_wikidata(x)
    }
  )
}
