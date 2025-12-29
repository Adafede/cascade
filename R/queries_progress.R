' Queries progress
#'
#' @description
#' Generates SPARQL queries for multiple taxa by fetching a query template
#' from a remote repository and parameterizing it with taxon QIDs and filters.
#'
#' @param xs Named list of taxon QIDs
#' @param start Start year for publication date filter (character)
#' @param end End year for publication date filter (character)
#' @param limit Maximum number of results per query (character)
#' @param query_url URL to the remote SPARQL query template. If NULL, uses default.
#' @param query_part_1 Deprecated. Kept for backward compatibility.
#' @param query_part_2 Deprecated. Kept for backward compatibility.
#' @param query_part_3 Deprecated. Kept for backward compatibility.
#' @param query_part_4 Deprecated. Kept for backward compatibility.
#'
#' @return A named list of parameterized SPARQL queries
#'
#' @examples
#' \dontrun{
#' qids <- list(Swertia = "Q1234", Kopsia = "Q5678")
#' queries <- queries_progress(xs = qids, start = "2000", end = "2024")
#' }
queries_progress <- function(
  xs,
  start = "0",
  end = "9999",
  limit = "1000000",
  query_url = NULL,
  query_part_1 = NULL,
  query_part_2 = NULL,
  query_part_3 = NULL,
  query_part_4 = NULL
) {
  # Backward compatibility: if old parameters are provided, use them
  if (
    !is.null(query_part_1) &&
      !is.null(query_part_2) &&
      !is.null(query_part_3) &&
      !is.null(query_part_4)
  ) {
    message(
      "Using legacy query construction (deprecated). Consider migrating to remote SPARQL queries."
    )

    return(
      xs |>
        purrr::map(
          .progress = TRUE,
          .f = function(
            x,
            start,
            end,
            limit,
            query_part_1,
            query_part_2,
            query_part_3,
            query_part_4
          ) {
            paste0(
              query_part_1,
              x,
              query_part_2,
              start,
              query_part_3,
              end,
              query_part_4,
              paste("\nLIMIT", limit)
            )
          },
          start = start,
          end = end,
          limit = limit,
          query_part_1 = query_part_1,
          query_part_2 = query_part_2,
          query_part_3 = query_part_3,
          query_part_4 = query_part_4
        )
    )
  }

  # Modern approach: fetch query from remote URL
  default_query_url <- "https://lotus.nprod.net/lotus-sparql-examples/examples/NPs/wd_nps_source_date_cascade_scholarly_subgraph.rq"
  query_url <- query_url %||% default_query_url

  # Fetch the query template once (efficiency)
  query_template <- tryCatch(
    {
      query_url |>
        readLines() |>
        paste0(collapse = "\n")
    },
    error = function(e) {
      stop(
        "Failed to fetch SPARQL query from: ",
        query_url,
        "\n",
        "Error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  # Generate queries for each taxon with progress tracking
  xs |>
    purrr::map(
      .progress = TRUE,
      .f = function(taxon_qid, template, start_year, end_year, result_limit) {
        # Substitute parameters in the template
        template |>
          gsub(
            pattern = "Q158572",
            replacement = taxon_qid,
            fixed = TRUE
          ) |>
          gsub(
            pattern = " 0 ",
            replacement = start_year,
            fixed = TRUE
          ) |>
          gsub(
            pattern = " 9999 ",
            replacement = end_year,
            fixed = TRUE
          ) |>
          gsub(
            pattern = "1000000",
            replacement = result_limit,
            fixed = TRUE
          )
      },
      template = query_template,
      start_year = start,
      end_year = end,
      result_limit = limit
    )
}
