#' Taxon name to QID
#'
#' @export
#'
#' @param taxon_name Taxon name
#'
#' @return A QID
#'
#' @examples
#' \dontrun{
#' taxon_name_to_qid(taxon_name = "Gentiana lutea")
#' }
taxon_name_to_qid <- function(taxon_name) {
  WikidataQueryServiceR::query_wikidata(
    sparql_query = paste0(
      "SELECT ?search ?item WHERE {
    SERVICE wikibase:mwapi {
      bd:serviceParam wikibase:endpoint \"www.wikidata.org\";
                      wikibase:api \"EntitySearch\";
                      mwapi:search \"",
      taxon_name,
      "                 \";
                      mwapi:language \"mul\".
    ?item wikibase:apiOutputItem mwapi:item.
    ?num wikibase:apiOrdinal true.
    }
    ?item (wdt:P225) ?search.
    FILTER (?num = 0)
    }"
    )
  ) |>
    tidytable::pull(item) |>
    gsub(pattern = "http://www.wikidata.org/entity/", replacement = "")
}
