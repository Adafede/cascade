#' Convert a scientific taxon name to its Wikidata QID
#'
#' Queries Wikidata for the item whose taxon name (`wdt:P225`) matches
#' `taxon_name` and returns the corresponding Wikidata entity identifier
#' (e.g., `"Q158572"`).
#'
#' @export
#'
#' @param taxon_name Character(1). Scientific name of the taxon (e.g.,
#'   `"Gentiana lutea"` or `"Swertia"`).
#'
#' @return Character vector of matching Wikidata QID(s).
#'
#' @concept EDAM:operation_0224 Query and retrieval
#' @concept EDAM:operation_3460 Taxonomic classification
#' @concept EDAM:data_1872 Taxonomic classification
#' @concept EDAM:topic_0637 Taxonomy
#'
#' @examples
#' \donttest{
#' try(taxon_name_to_qid(taxon_name = "Gentiana lutea"))
#' }
taxon_name_to_qid <- function(taxon_name) {
  valid <- is.character(taxon_name) &&
    length(taxon_name) == 1L &&
    !is.na(taxon_name)
  stopifnot(
    "`taxon_name` must be a single non-missing character string." = valid
  )

  query <- sprintf(
    paste0(
      "PREFIX wdt: <http://www.wikidata.org/prop/direct/>\n",
      "SELECT ?taxon WHERE {\n",
      "  ?taxon wdt:P225 \"%s\" .\n",
      "}"
    ),
    gsub('"', '\\"', taxon_name, fixed = TRUE)
  )

  res <- query_wikidata(sparql_query = query)
  if (!("taxon" %in% names(res)) || nrow(res) == 0L) {
    return(character(0L))
  }
  sub("^https?://www\\.wikidata\\.org/entity/", "", res$taxon)
}
