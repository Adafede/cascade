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
  # Fetch the SPARQL query template from remote repository
  query_template <- "https://adafede.github.io/sparql-examples/examples/Taxa/wd_taxa_name_to_qid.rq" |>
    readLines() |>
    paste0(collapse = "\n")

  # Substitute the taxon name parameter
  query <- gsub(
    pattern = "Gentiana lutea",
    replacement = taxon_name,
    x = query_template,
    fixed = FALSE
  )

  # Execute query and extract QID
  query_wikidata(sparql_query = query) |>
    tidytable::pull(taxon) |>
    gsub(pattern = "http://www.wikidata.org/entity/", replacement = "")
}
