if (!requireNamespace("httr2", quietly = TRUE)) {
  install.packages("httr2")
}
if (!requireNamespace("tidytable", quietly = TRUE)) {
  install.packages("tidytable")
}

#' @title Query WDQS (Wikidata Query Service)
#'
#' @description
#' Performs a SPARQL query against the Wikidata Query Service (WDQS) or another compatible endpoint, returning results as a tidytable. Handles datetime conversion and optional URL prefix removal.
#'
#' @export
#' 
#' @param sparql_query Character. SPARQL query string.
#' @param remove_url Logical. If TRUE, removes 'http://www.wikidata.org/entity/' prefix from character columns.
#' @param endpoint Character. SPARQL endpoint URL.
#' @param agent Character. User agent string for polite querying.
#' @param headers Character. HTTP Accept header for the request.
#' @param fallback Logical. If TRUE, retries with qlever.dev endpoint on failure.
#'
#' @return A tidytable::tidytable data frame containing the results of the query.
#'
#' @examples
#' NULL
query_wikidata <- function(
    sparql_query,
    remove_url = TRUE,
    endpoint = "https://query.wikidata.org/sparql",
    agent = "https://github.com/bearloga/WikidataQueryServiceR",
    headers = "application/sparql-results+json",
    fallback = TRUE
) {
  # Input validation
  if (!is.character(sparql_query) || length(sparql_query) != 1) {
    stop("query must be a single character string (SPARQL query)")
  }
  if (!is.logical(remove_url) || length(remove_url) != 1) {
    stop("remove_url must be a single logical value")
  }
  if (!is.character(endpoint) || length(endpoint) != 1) {
    stop("endpoint must be a single character string (URL)")
  }
  if (!is.character(agent) || length(agent) != 1) {
    stop("agent must be a single character string")
  }
  if (!is.character(headers) || length(headers) != 1) {
    stop("headers must be a single character string")
  }
  if (!is.logical(fallback) || length(fallback) != 1) {
    stop("fallback must be a single logical value")
  }
  
  # Internal function to execute the query
  execute_query <- function(ep) {
    temp <- httr2::request(ep) |>
      httr2::req_url_query(query = sparql_query) |>
      httr2::req_user_agent(string = agent) |>
      httr2::req_headers(Accept = headers) |>
      httr2::req_perform() |>
      httr2::resp_body_json()
    
    # Handle empty result set
    if (length(temp$results$bindings) == 0) {
      return(tidytable::tidytable(matrix(
        character(),
        nrow = 0,
        ncol = length(temp$head$vars),
        dimnames = list(c(), unlist(temp$head$vars))
      )))
    }
    
    # Parse bindings into a character matrix
    bindings <- temp$results$bindings
    n_rows <- length(bindings)
    col_names <- temp$head$vars
    n_cols <- length(col_names)
    result_matrix <- matrix(
      character(n_rows * n_cols),
      nrow = n_rows,
      ncol = n_cols,
      dimnames = list(NULL, col_names)
    )
    for (j in seq_len(n_cols)) {
      col_name <- col_names[[j]]
      for (i in seq_len(n_rows)) {
        binding <- bindings[[i]][[col_name]]
        if (!is.null(binding)) {
          result_matrix[i, j] <- binding$value
        }
      }
    }
    df <- result_matrix |>
      tidytable::as_tidytable()
    
    # Convert datetime columns if present
    if (n_rows > 0) {
      first_binding <- bindings[[1]]
      datetime_cols <- character(0)
      for (col_name in col_names) {
        binding <- first_binding[[col_name]]
        if (
          !is.null(binding) &&
          !is.null(binding$datatype) &&
          binding$datatype == "http://www.w3.org/2001/XMLSchema#dateTime"
        ) {
          datetime_cols <- c(datetime_cols, col_name)
        }
      }
      if (length(datetime_cols) > 0) {
        df <- df |>
          tidytable::mutate(
            tidytable::across(
              tidytable::all_of(datetime_cols),
              .fns = as.POSIXct,
              format = "%Y-%m-%dT%H:%M:%SZ",
              tz = "GMT"
            )
          )
      }
    }
    
    # Optionally remove Wikidata entity URL prefix
    if (remove_url) {
      char_cols <- sapply(df, is.character)
      if (any(char_cols)) {
        df <- df |>
          tidytable::mutate(
            tidytable::across(
              tidytable::where(is.character),
              .fns = stringi::stri_replace_first_regex,
              pattern = "^http://www\\.wikidata\\.org/entity/",
              replacement = ""
            )
          )
      }
    }
    
    return(df)
  }
  
  # Try primary endpoint, with fallback if enabled
  result <- tryCatch(
    {
      execute_query(endpoint)
    },
    error = function(e) {
      if (fallback && endpoint != "https://qlever.dev/api/wikidata") {
        warning(
          "Primary endpoint failed: ",
          e$message,
          "\nRetrying with fallback endpoint: https://qlever.dev/api/wikidata",
          call. = FALSE
        )
        tryCatch(
          {
            execute_query("https://qlever.dev/api/wikidata")
          },
          error = function(e2) {
            stop(
              "Both primary and fallback endpoints failed.\n",
              "Primary error: ",
              e$message,
              "\n",
              "Fallback error: ",
              e2$message,
              call. = FALSE
            )
          }
        )
      } else {
        stop(e$message, call. = FALSE)
      }
    }
  )
  
  return(result)
}
