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
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      WikidataQueryServiceR::query_wikidata(x)
    }
  )
}