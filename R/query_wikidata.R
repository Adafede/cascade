if (!requireNamespace("curl", quietly = TRUE)) {
  install.packages("curl")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

#' @title Query a SPARQL endpoint efficiently
#'
#' @description
#' Performs a SPARQL query and returns a data.table. Optimised for very large
#' result sets (millions of rows):
#'
#'  - Requests **CSV** (no IRI angle-bracket decoration, smallest text format)
#'  - Requests **gzip** transfer encoding (3-5x less data over the wire)
#'  - Streams the response **directly to disk** via curl (zero R memory use
#'    during download)
#'  - Parses with **data.table::fread** (C-level, multi-threaded)
#'
#' Falls back to JSON for endpoints that do not support CSV.
#'
#' @param sparql_query  Character. SPARQL query string.
#' @param remove_url    Logical. Strip `http://www.wikidata.org/entity/` prefix
#'                      from character columns (default TRUE).
#' @param endpoint      Character. SPARQL endpoint URL.
#' @param agent         Character. User-Agent header string.
#' @param timeout       Integer. Total request timeout in seconds (default 3600).
#' @param fallback      Logical. Retry with QLever Wikidata endpoint on failure.
#' @param headers       Character or NULL. Optional `Accept` header used for
#'                      backward compatibility with older versions. If set,
#'                      this value is used as preferred response format.
#'
#' @return A `data.table`.
#' @export
#' @examples NULL
query_wikidata <- function(
  sparql_query,
  remove_url = TRUE,
  endpoint = "https://query.wikidata.org/sparql",
  agent = "https://github.com/bearloga/WikidataQueryServiceR",
  timeout = 3600L,
  fallback = TRUE,
  headers = NULL
) {
  # Support legacy signature where `headers` was argument #5.
  if (is.character(timeout) && length(timeout) == 1L && is.null(headers)) {
    headers <- timeout
    timeout <- 3600L
  }

  # Support accidental positional use with both legacy and current conventions:
  # query_wikidata(..., <timeout>, <fallback>) while still allowing `headers`.
  if (
    !is.null(headers) &&
      (is.numeric(headers)) &&
      (is.logical(timeout) || is.numeric(timeout) || is.integer(timeout))
  ) {
    fallback <- as.logical(timeout)[1L]
    timeout <- headers
    headers <- NULL
  }

  # ── Validation ──────────────────────────────────────────────────────────────
  if (
    !is.character(sparql_query) ||
      length(sparql_query) != 1L ||
      is.na(sparql_query)
  ) {
    stop(
      "`sparql_query` must be a single non-missing character string.",
      call. = FALSE
    )
  }
  if (
    !is.logical(remove_url) || length(remove_url) != 1L || is.na(remove_url)
  ) {
    stop(
      "`remove_url` must be a single non-missing logical value.",
      call. = FALSE
    )
  }
  if (!is.character(endpoint) || length(endpoint) != 1L || is.na(endpoint)) {
    stop(
      "`endpoint` must be a single non-missing character string.",
      call. = FALSE
    )
  }
  if (!is.character(agent) || length(agent) != 1L || is.na(agent)) {
    stop(
      "`agent` must be a single non-missing character string.",
      call. = FALSE
    )
  }
  if (
    !(is.null(headers) ||
      (is.character(headers) && length(headers) == 1L && !is.na(headers)))
  ) {
    stop(
      "`headers` must be NULL or a single non-missing character string.",
      call. = FALSE
    )
  }
  if (
    !(is.numeric(timeout)) ||
      length(timeout) != 1L ||
      is.na(timeout) ||
      timeout <= 0
  ) {
    stop(
      "`timeout` must be a single positive numeric value (seconds).",
      call. = FALSE
    )
  }
  if (!is.logical(fallback) || length(fallback) != 1L || is.na(fallback)) {
    stop(
      "`fallback` must be a single non-missing logical value.",
      call. = FALSE
    )
  }

  # ── Helpers ─────────────────────────────────────────────────────────────────

  # Build the full request URL with query param
  build_url <- function(ep) {
    sep <- if (grepl("?", ep, fixed = TRUE)) "&" else "?"
    paste0(ep, sep, "query=", curl::curl_escape(sparql_query))
  }

  # Strip Wikidata entity prefix from all character columns (vectorised)
  strip_wd <- function(dt) {
    wd <- "^https?://www\\.wikidata\\.org/entity/"
    char_j <- which(vapply(dt, is.character, logical(1L)))
    for (j in char_j) {
      v <- dt[[j]]
      if (any(grepl(wd, v, perl = TRUE), na.rm = TRUE)) {
        data.table::set(dt, j = j, value = sub(wd, "", v, perl = TRUE))
      }
    }
    dt
  }

  # ── Primary path: CSV + gzip via curl ───────────────────────────────────────
  execute_csv <- function(ep, accept = "text/csv") {
    tmp <- tempfile(fileext = ".csv")
    on.exit(unlink(tmp), add = TRUE)

    h <- curl::new_handle()
    curl::handle_setheaders(
      h,
      "Accept" = accept,
      "Accept-Encoding" = "gzip, deflate", # compress transfer
      "User-Agent" = agent
    )
    curl::handle_setopt(h, timeout = as.integer(timeout), followlocation = 1L)

    # Stream to disk — curl decompresses gzip transparently
    curl::curl_download(build_url(ep), destfile = tmp, handle = h, quiet = TRUE)

    # Detect if the server returned an error document instead of CSV
    # (some endpoints return 200 + HTML error page)
    first_line <- readLines(tmp, n = 1L, warn = FALSE)
    if (length(first_line) && startsWith(trimws(first_line), "<")) {
      stop(
        "Endpoint returned non-CSV response (HTML/XML). ",
        "Check the URL or query syntax.",
        call. = FALSE
      )
    }

    # fread: C-level, multi-threaded, handles gzip natively too
    dt <- data.table::fread(
      tmp,
      sep = ",",
      header = TRUE,
      encoding = "UTF-8",
      showProgress = FALSE,
      data.table = TRUE
    )

    if (remove_url && nrow(dt) > 0L) {
      dt <- strip_wd(dt)
    }
    dt
  }

  # ── JSON fallback (endpoints without CSV support) ───────────────────────────
  execute_json <- function(ep, accept = "application/sparql-results+json") {
    if (!requireNamespace("httr2", quietly = TRUE)) {
      install.packages("httr2")
    }

    resp <- httr2::request(ep) |>
      httr2::req_url_query(query = sparql_query) |>
      httr2::req_user_agent(agent) |>
      httr2::req_headers(
        Accept = accept,
        `Accept-Encoding` = "gzip, deflate"
      ) |>
      httr2::req_timeout(timeout) |>
      httr2::req_perform()

    body <- httr2::resp_body_json(resp)
    col_names <- unlist(body$head$vars)
    bindings <- body$results$bindings

    if (length(bindings) == 0L) {
      return(data.table::as.data.table(
        stats::setNames(
          replicate(length(col_names), character(0L), simplify = FALSE),
          col_names
        )
      ))
    }

    # Vectorised column extraction — no nested for-loops
    dt <- data.table::as.data.table(
      stats::setNames(
        lapply(col_names, \(col) {
          vapply(
            bindings,
            \(row) {
              v <- row[[col]]$value
              if (is.null(v)) NA_character_ else v
            },
            character(1L)
          )
        }),
        col_names
      )
    )

    # xsd:dateTime detection & conversion
    for (col in col_names) {
      is_datetime <- any(vapply(
        bindings,
        FUN = function(row) {
          dt_type <- row[[col]]$datatype
          !is.null(dt_type) &&
            dt_type == "http://www.w3.org/2001/XMLSchema#dateTime"
        },
        FUN.VALUE = logical(1L)
      ))
      if (is_datetime) {
        data.table::set(
          dt,
          j = col,
          value = as.POSIXct(
            dt[[col]],
            format = "%Y-%m-%dT%H:%M:%SZ",
            tz = "GMT"
          )
        )
      }
    }

    if (remove_url && nrow(dt) > 0L) {
      dt <- strip_wd(dt)
    }
    dt
  }

  # ── Execute with optional fallback ──────────────────────────────────────────
  execute <- function(ep) {
    prefer_json <- !is.null(headers) &&
      grepl("json", headers, ignore.case = TRUE)
    preferred <- if (prefer_json) execute_json else execute_csv
    alternate <- if (identical(preferred, execute_csv)) {
      execute_json
    } else {
      execute_csv
    }

    preferred_accept <- headers %||% if (identical(preferred, execute_csv)) {
        "text/csv"
      } else {
        "application/sparql-results+json"
      }
    alternate_accept <- if (identical(alternate, execute_csv)) {
      "text/csv"
    } else {
      "application/sparql-results+json"
    }

    tryCatch(
      preferred(ep, accept = preferred_accept),
      error = function(e) {
        message(
          "Preferred format attempt failed on ",
          ep,
          " (",
          conditionMessage(e),
          ") - retrying with alternate format."
        )
        alternate(ep, accept = alternate_accept)
      }
    )
  }

  fallback_ep <- "https://qlever.dev/api/wikidata"

  tryCatch(
    execute(endpoint),
    error = function(e) {
      if (fallback && !identical(endpoint, fallback_ep)) {
        warning(
          "Primary endpoint failed: ",
          conditionMessage(e),
          "\nRetrying with fallback: ",
          fallback_ep,
          call. = FALSE
        )
        tryCatch(
          execute(fallback_ep),
          error = function(e2) {
            stop(
              "Both endpoints failed.\n",
              "Primary : ",
              conditionMessage(e),
              "\n",
              "Fallback: ",
              conditionMessage(e2),
              call. = FALSE
            )
          }
        )
      } else {
        stop(conditionMessage(e), call. = FALSE)
      }
    }
  )
}
