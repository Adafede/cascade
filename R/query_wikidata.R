#' Query a SPARQL endpoint efficiently
#'
#' Performs a SPARQL query against Wikidata (or another SPARQL 1.1 endpoint)
#' and returns a `data.table`. Requests gzip-compressed CSV (falling back to
#' TSV) streamed directly to disk via `curl` and parsed with
#' [data.table::fread()]. Automatically retries against the QLever Wikidata
#' endpoint (`https://qlever.dev/api/wikidata`) when `fallback = TRUE`.
#'
#' @param sparql_query Character(1). SPARQL query string.
#' @param remove_url Logical(1). Strip `http(s)://www.wikidata.org/entity/`
#'   prefix from character columns (default `TRUE`).
#' @param endpoint Character(1). SPARQL endpoint URL. Default is
#'   `"https://query.wikidata.org/sparql"`.
#' @param agent Character(1). HTTP `User-Agent` header string.
#' @param timeout Numeric(1). Total request timeout in seconds (default
#'   `3600L`).
#' @param fallback Logical(1). Retry with QLever Wikidata endpoint on failure
#'   (default `TRUE`).
#' @param headers Character(1) or `NULL`. Optional `Accept` header for backward
#'   compatibility.
#' @param post Logical(1). Send the SPARQL query as an HTTP POST request body
#'   (`application/x-www-form-urlencoded`) rather than a GET query parameter
#'   (default `FALSE`).
#'
#' @return A `data.table` containing the query results.
#' @export
#'
#' @concept EDAM:operation_0224 Query and retrieval
#' @concept EDAM:data_1086 Compound identifier
#' @concept EDAM:format_3790 SPARQL
#' @concept EDAM:topic_3172 Metabolomics
#'
#' @examples
#' \donttest{
#' q <- paste0(
#'   "SELECT ?taxon WHERE { ",
#'   "?taxon <http://www.wikidata.org/prop/direct/P225> 'Gentiana lutea' . }"
#' )
#' try(query_wikidata(q, timeout = 15L))
#' }
query_wikidata <- function(
  sparql_query,
  remove_url = TRUE,
  endpoint = "https://query.wikidata.org/sparql",
  agent = "https://github.com/adafede/cascade",
  timeout = 3600L,
  fallback = TRUE,
  headers = NULL,
  post = FALSE
) {
  if (is.character(timeout) && length(timeout) == 1L && is.null(headers)) {
    headers <- timeout
    timeout <- 3600L
  }

  if (
    !is.null(headers) &&
      is.numeric(headers) &&
      (is.logical(timeout) || is.numeric(timeout) || is.integer(timeout))
  ) {
    fallback <- as.logical(timeout)[1L]
    timeout <- headers
    headers <- NULL
  }

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
  valid_headers <- is.null(headers) ||
    (is.character(headers) && length(headers) == 1L && !is.na(headers))
  if (!valid_headers) {
    stop(
      "`headers` must be NULL or a single non-missing character string.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(timeout) ||
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
  if (!is.logical(post) || length(post) != 1L || is.na(post)) {
    stop(
      "`post` must be a single non-missing logical value.",
      call. = FALSE
    )
  }

  build_url <- function(ep) {
    sep <- if (grepl("?", ep, fixed = TRUE)) "&" else "?"
    paste0(ep, sep, "query=", curl::curl_escape(sparql_query))
  }

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

  execute_delimited <- function(
    ep,
    accept = "text/csv",
    sep = ",",
    use_post = post
  ) {
    tmp <- tempfile(fileext = if (sep == "\t") ".tsv" else ".csv")
    on.exit(unlink(tmp), add = TRUE)

    h <- curl::new_handle()
    curl::handle_setheaders(
      h,
      "Accept" = accept,
      "Accept-Encoding" = "gzip, deflate",
      "User-Agent" = agent
    )
    curl::handle_setopt(h, timeout = as.integer(timeout), followlocation = 1L)

    if (use_post) {
      curl::handle_setopt(
        h,
        post = 1L,
        postfields = paste0("query=", curl::curl_escape(sparql_query))
      )
      curl::handle_setheaders(
        h,
        "Content-Type" = "application/x-www-form-urlencoded"
      )
      url <- ep
    } else {
      url <- build_url(ep)
    }

    curl::curl_download(url, destfile = tmp, handle = h, quiet = TRUE)

    first_line <- readLines(tmp, n = 1L, warn = FALSE)
    if (length(first_line) && startsWith(trimws(first_line), "<")) {
      stop(
        "Endpoint returned non-tabular response (HTML/XML). ",
        "Check the URL or query syntax.",
        call. = FALSE
      )
    }

    dt <- data.table::fread(
      tmp,
      sep = sep,
      quote = if (sep == "\t") "" else "\"",
      header = TRUE,
      encoding = "UTF-8",
      showProgress = FALSE,
      data.table = TRUE
    )

    if (sep == "\t") {
      data.table::setnames(dt, sub("^\\?", "", names(dt)))
      if (nrow(dt) > 0L) {
        for (j in which(vapply(dt, is.character, logical(1L)))) {
          v <- dt[[j]]
          v <- sub("^<(.+)>$", "\\1", v)
          v <- sub('^"(.*)"(@[a-zA-Z0-9-]+|\\^\\^<[^>]+>)?$', "\\1", v)
          v <- gsub('\\"', '"', v, fixed = TRUE)
          data.table::set(dt, j = j, value = v)
        }
      }
    }

    if (remove_url && nrow(dt) > 0L) {
      dt <- strip_wd(dt)
    }
    dt
  }

  execute <- function(ep, use_post = post) {
    tryCatch(
      execute_delimited(
        ep,
        accept = "text/csv",
        sep = ",",
        use_post = use_post
      ),
      error = function(e) {
        message(
          "CSV attempt failed on ",
          ep,
          " (",
          conditionMessage(e),
          ") - retrying with TSV format."
        )
        execute_delimited(
          ep,
          accept = "text/tab-separated-values",
          sep = "\t",
          use_post = use_post
        )
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
          execute(fallback_ep, use_post = TRUE),
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

#' Build parameterized SPARQL queries for a list of taxon QIDs
#'
#' @param xs Named list of taxon QIDs.
#' @param start Start year for publication date filter (character).
#' @param end End year for publication date filter (character).
#' @param limit Maximum number of results per query (character).
#' @param query_url Optional URL to a remote SPARQL query template.
#' @param query_part_1 Optional legacy query prefix.
#' @param query_part_2 Optional legacy query segment 2.
#' @param query_part_3 Optional legacy query segment 3.
#' @param query_part_4 Optional legacy query suffix.
#'
#' @return A named list of SPARQL query strings.
#' @keywords internal
#' @noRd
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
  if (
    !is.null(query_part_1) &&
      !is.null(query_part_2) &&
      !is.null(query_part_3) &&
      !is.null(query_part_4)
  ) {
    return(lapply(xs, function(x) {
      paste0(
        query_part_1,
        x,
        query_part_2,
        start,
        query_part_3,
        end,
        query_part_4,
        "\nLIMIT ",
        limit
      )
    }))
  }

  if (!is.null(query_url)) {
    query_template <- tryCatch(
      paste0(
        suppressWarnings(readLines(query_url, warn = FALSE)),
        collapse = "\n"
      ),
      error = function(e) {
        stop(
          "Failed to fetch SPARQL query from: ",
          query_url,
          "\nError: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  } else {
    query_template <- paste0(
      "SELECT DISTINCT ?structure ?structureLabel ?structure_id ",
      "?structureSmiles ?taxaLabels ?taxa ?referencesLabels ",
      "?references_ids ?references WHERE {\n",
      "  ?taxa (wdt:P171*) wd:Q158572;\n",
      "  wdt:P225 ?taxaLabels.\n",
      "  ?structure wdt:P235 ?structure_id;\n",
      "  wdt:P233 ?structureSmiles;\n",
      "  p:P703 ?statement.\n",
      "  ?statement ps:P703 ?taxa;\n",
      "  prov:wasDerivedFrom ?ref.\n",
      "  ?ref pr:P248 ?references.\n",
      "  SERVICE <https://query-scholarly.wikidata.org/sparql> {\n",
      "    ?references wdt:P1476 ?referencesLabels;\n",
      "    wdt:P356 ?references_ids;\n",
      "    wdt:P577 ?art_date.\n",
      "  }\n",
      "  FILTER(((YEAR(?art_date)) >= 0 ) && ((YEAR(?art_date)) <= 9999 ))\n",
      "  SERVICE wikibase:label { ",
      "bd:serviceParam wikibase:language \"[AUTO_LANGUAGE], mul, en\". }\n",
      "}\nLIMIT 1000000"
    )
  }

  lapply(xs, function(taxon_qid) {
    q <- gsub("Q158572", taxon_qid, query_template, fixed = TRUE)
    q <- gsub(" 0 ", paste0(" ", start, " "), q, fixed = TRUE)
    q <- gsub(" 9999 ", paste0(" ", end, " "), q, fixed = TRUE)
    gsub("1000000", limit, q, fixed = TRUE)
  })
}

#' Execute a list of SPARQL queries against Wikidata
#'
#' @param xs Named list of SPARQL query strings.
#'
#' @return A named list of `data.table` query results.
#' @keywords internal
#' @noRd
wiki_progress <- function(xs) {
  lapply(xs, query_wikidata)
}

#' Join Wikidata query results with LOTUS chemical classifications
#'
#' @param xs Named list of `data.table` results from `wiki_progress()`.
#' @param structures_classified `data.table` of LOTUS classified structures.
#'
#' @return A named list of cleaned and classified `data.table` objects.
#' @keywords internal
#' @noRd
tables_progress <- function(xs, structures_classified) {
  sc_dt <- data.table::as.data.table(structures_classified)

  lapply(xs, function(x) {
    dt <- data.table::copy(data.table::as.data.table(x))
    if (nrow(dt) > 0L) {
      joined <- dt_left_join(dt, sc_dt)
      joined[,
        structureImage := encode_smiles(
          structureSmiles,
          empty = NA_character_,
          reserved = TRUE
        )
      ]

      drop_cols <- intersect(
        c("references_ids", "structure_id", "structureSmiles"),
        names(joined)
      )
      if (length(drop_cols) > 0L) {
        joined[, (drop_cols) := NULL]
      }

      first_cols <- intersect(
        c("structureLabel", "structure", "structureImage"),
        names(joined)
      )
      data.table::setcolorder(
        joined,
        c(first_cols, setdiff(names(joined), first_cols))
      )

      for (grp_col in c(
        "structureLabel",
        "chemical_class",
        "chemical_superclass",
        "chemical_pathway"
      )) {
        if (grp_col %in% names(joined)) {
          joined[, .tmp_n := .N, by = grp_col]
          data.table::setorder(joined, -.tmp_n, na.last = TRUE)
          joined[, .tmp_n := NULL]
        }
      }

      unique(joined)
    } else {
      data.table::data.table(
        structureLabel = NA_character_,
        structure = NA_character_,
        structureImage = NA_character_,
        taxaLabels = NA_character_,
        taxa = NA_character_,
        referencesLabels = NA_character_,
        references = NA_character_,
        chemical_pathway = NA_character_,
        chemical_superclass = NA_character_,
        chemical_class = NA_character_
      )
    }
  })
}

#' Format a Wikidata compound-taxon-reference table as an interactive gt table
#'
#' Converts `structure`, `taxa`, and `references` columns into clickable HTML
#' links and renders `structureImage` SMILES strings as 2D chemical depictions
#' via `molinfo()`.
#'
#' @param table Data frame or `data.table` from `tables_progress()`.
#' @param title Table header title (default `""`).
#' @param subtitle Table header subtitle (default `""`).
#'
#' @return A `gt_tbl` object.
#' @keywords internal
#' @noRd
format_gt <- function(table, title = "", subtitle = "") {
  df <- as.data.frame(table, stringsAsFactors = FALSE)
  make_wd_links <- function(vals) {
    lapply(as.character(vals), function(v) {
      if (is.na(v) || !nzchar(v)) {
        return(gt::html(""))
      }
      href <- if (grepl("^https?://", v)) {
        v
      } else {
        paste0("http://www.wikidata.org/entity/", v)
      }
      gt::html(sprintf('<a href="%s">%s</a>', href, v))
    })
  }

  for (col in intersect(c("structure", "taxa", "references"), names(df))) {
    df[[col]] <- I(make_wd_links(df[[col]]))
  }

  gt_out <- gt::gt(df) |>
    gt::tab_header(
      title = gt::md(title),
      subtitle = gt::md(subtitle)
    )

  if ("structureImage" %in% names(df)) {
    gt_out <- gt::text_transform(
      gt_out,
      locations = gt::cells_body(columns = structureImage),
      fn = molinfo
    )
  }

  gt_out
}
