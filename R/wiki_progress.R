#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
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
