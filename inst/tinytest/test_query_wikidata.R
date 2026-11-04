library(tinytest)
library(cascade)

queries_progress <- utils::getFromNamespace("queries_progress", "cascade")
wiki_progress <- utils::getFromNamespace("wiki_progress", "cascade")
tables_progress <- utils::getFromNamespace("tables_progress", "cascade")
format_gt <- utils::getFromNamespace("format_gt", "cascade")

# --- 1. query_wikidata input validation ---
expect_error(
  query_wikidata(character(0L)),
  pattern = "`sparql_query` must be a single non-missing character string"
)
expect_error(
  query_wikidata(NA_character_),
  pattern = "`sparql_query` must be a single non-missing character string"
)
expect_error(
  query_wikidata("SELECT * WHERE {}", remove_url = NA),
  pattern = "`remove_url` must be a single non-missing logical value"
)
expect_error(
  query_wikidata("SELECT * WHERE {}", endpoint = 123),
  pattern = "`endpoint` must be a single non-missing character string"
)
expect_error(
  query_wikidata("SELECT * WHERE {}", agent = NA_character_),
  pattern = "`agent` must be a single non-missing character string"
)
expect_error(
  query_wikidata("SELECT * WHERE {}", headers = c("a", "b")),
  pattern = "`headers` must be NULL or a single non-missing character string"
)
expect_error(
  query_wikidata("SELECT * WHERE {}", timeout = 0),
  pattern = "`timeout` must be a single positive numeric value"
)
expect_error(
  query_wikidata("SELECT * WHERE {}", fallback = NA),
  pattern = "`fallback` must be a single non-missing logical value"
)
expect_error(
  query_wikidata("SELECT * WHERE {}", post = "TRUE"),
  pattern = "`post` must be a single non-missing logical value"
)

# --- 2. query_wikidata offline execution & branches via mocked curl_download ---
local({
  curl_ns <- asNamespace("curl")
  orig_dl <- get("curl_download", envir = curl_ns)
  unlockBinding("curl_download", curl_ns)
  on.exit(
    {
      assign("curl_download", orig_dl, envir = curl_ns)
      lockBinding("curl_download", curl_ns)
    },
    add = TRUE
  )

# 2a. CSV success with remove_url = TRUE vs remove_url = FALSE, GET with '?' in endpoint
last_url <- NULL
assign(
  "curl_download",
  function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
    last_url <<- url
    writeLines(
      c(
        "taxon,label,count",
        "http://www.wikidata.org/entity/Q158572,Swertia,10",
        "https://www.wikidata.org/entity/Q12345,Gentiana,20"
      ),
      destfile
    )
    invisible(destfile)
  },
  envir = curl_ns
)

res_strip <- query_wikidata(
  "SELECT * WHERE {}",
  remove_url = TRUE,
  endpoint = "https://example.org/sparql?format=csv",
  fallback = FALSE
)
expect_true(grepl("&query=", last_url, fixed = TRUE))
expect_identical(res_strip$taxon, c("Q158572", "Q12345"))
expect_identical(res_strip$label, c("Swertia", "Gentiana"))
expect_identical(res_strip$count, c(10L, 20L))

res_keep_url <- query_wikidata(
  "SELECT * WHERE {}",
  remove_url = FALSE,
  endpoint = "https://example.org/sparql",
  post = TRUE,
  fallback = FALSE
)
expect_identical(last_url, "https://example.org/sparql")
expect_identical(
  res_keep_url$taxon,
  c(
    "http://www.wikidata.org/entity/Q158572",
    "https://www.wikidata.org/entity/Q12345"
  )
)

# 2b. Legacy argument shims: character timeout -> headers, numeric headers -> timeout
res_shim1 <- query_wikidata(
  "SELECT * WHERE {}",
  timeout = "text/csv",
  fallback = FALSE
)
expect_identical(nrow(res_shim1), 2L)

res_shim2 <- query_wikidata(
  "SELECT * WHERE {}",
  timeout = FALSE,
  headers = 30L
)
expect_identical(nrow(res_shim2), 2L)

# 2c. HTML/XML error response detection + TSV fallback with RDF literal/URI stripping
call_idx <- 0L
assign(
  "curl_download",
  function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
    call_idx <<- call_idx + 1L
    if (call_idx == 1L) {
      writeLines("  <html><body>Error 502</body></html>", destfile)
    } else {
      writeLines(
        c(
          "?taxon\t?label",
          "<http://www.wikidata.org/entity/Q158572>\tSwertia",
          "<https://www.wikidata.org/entity/Q999>\tGentiana"
        ),
        destfile
      )
    }
    invisible(destfile)
  },
  envir = curl_ns
)
res_tsv_fb <- query_wikidata("SELECT ?taxon ?label WHERE {}", fallback = FALSE)
expect_identical(names(res_tsv_fb), c("taxon", "label"))
expect_identical(res_tsv_fb$taxon, c("Q158572", "Q999"))
expect_identical(res_tsv_fb$label, c("Swertia", "Gentiana"))

# 2d. Primary endpoint fails -> warns and succeeds on QLever fallback endpoint
assign(
  "curl_download",
  function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
    if (!grepl("qlever", url, fixed = TRUE)) {
      stop("Primary connection timeout")
    }
    writeLines(c("taxon", "http://www.wikidata.org/entity/Q777"), destfile)
    invisible(destfile)
  },
  envir = curl_ns
)
expect_warning(
  res_qlever <- query_wikidata("SELECT ?taxon WHERE {}", fallback = TRUE),
  pattern = "Primary endpoint failed.*Retrying with fallback"
)
expect_identical(res_qlever$taxon, "Q777")

# 2e. Both primary and fallback endpoints fail -> stops with combined message
assign(
  "curl_download",
  function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
    stop("Network unreachable")
  },
  envir = curl_ns
)
expect_warning(
  expect_error(
    query_wikidata("SELECT ?taxon WHERE {}", fallback = TRUE),
    pattern = "Both endpoints failed"
  ),
  pattern = "Primary endpoint failed"
)
expect_error(
  query_wikidata("SELECT ?taxon WHERE {}", fallback = FALSE),
  pattern = "Network unreachable"
)

# --- 3. taxon_name_to_qid ---
expect_error(
  taxon_name_to_qid(123),
  pattern = "`taxon_name` must be a single non-missing character string"
)
expect_error(
  taxon_name_to_qid(c("Swertia", "Gentiana")),
  pattern = "`taxon_name` must be a single non-missing character string"
)
expect_error(
  taxon_name_to_qid(NA_character_),
  pattern = "`taxon_name` must be a single non-missing character string"
)

# Successful QID resolution and empty result handling
assign(
  "curl_download",
  function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
    writeLines(c("taxon", "http://www.wikidata.org/entity/Q158572"), destfile)
    invisible(destfile)
  },
  envir = curl_ns
)
expect_identical(taxon_name_to_qid("Swertia"), "Q158572")

assign(
  "curl_download",
  function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
    writeLines("taxon", destfile)
    invisible(destfile)
  },
  envir = curl_ns
)
expect_identical(taxon_name_to_qid("NonExistentTaxon"), character(0L))

# --- 4. wiki_progress ---
assign(
  "curl_download",
  function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
    writeLines(c("structure", "Q415869"), destfile)
    invisible(destfile)
  },
  envir = curl_ns
)
  wp_out <- wiki_progress(list(Swertia = "SELECT ?structure WHERE {}"))
  expect_identical(names(wp_out), "Swertia")
  expect_identical(wp_out$Swertia$structure, "Q415869")
})

# --- 5. queries_progress (default template, custom query_url, failing query_url, legacy parts) ---
local({
  # 5a. Default template
  qp_def <- queries_progress(
    list(Swertia = "Q158572", Gentiana = "Q158500"),
    start = "1950",
    end = "2025",
    limit = "500"
  )
  expect_identical(names(qp_def), c("Swertia", "Gentiana"))
  expect_true(grepl("wd:Q158572;", qp_def$Swertia, fixed = TRUE))
  expect_true(grepl("wd:Q158500;", qp_def$Gentiana, fixed = TRUE))
  expect_true(grepl(">= 1950 )", qp_def$Swertia, fixed = TRUE))
  expect_true(grepl("<= 2025 )", qp_def$Swertia, fixed = TRUE))
  expect_true(grepl("LIMIT 500$", qp_def$Swertia))

  # 5b. Custom query_url from local file
  tmpl_file <- tempfile(fileext = ".sparql")
  on.exit(unlink(tmpl_file), add = TRUE)
  writeLines(
    c(
      "SELECT * WHERE { ?x wdt:P171* wd:Q158572 . FILTER(YEAR(?d) >= 0 && YEAR(?d) <= 9999 ) }",
      "LIMIT 1000000"
    ),
    tmpl_file
  )
  qp_url <- queries_progress(
    list(Test = "Q99999"),
    start = "2000",
    end = "2020",
    limit = "10",
    query_url = tmpl_file
  )
  expect_true(grepl("wd:Q99999", qp_url$Test, fixed = TRUE))
  expect_true(grepl(">= 2000 ", qp_url$Test, fixed = TRUE))
  expect_true(grepl("<= 2020 ", qp_url$Test, fixed = TRUE))
  expect_true(grepl("LIMIT 10$", qp_url$Test))

  # 5c. Failing query_url throws informative error
  expect_error(
    queries_progress(
      list(Test = "Q1"),
      query_url = "/nonexistent/path/to/query.sparql"
    ),
    pattern = "Failed to fetch SPARQL query from"
  )

  # 5d. Legacy 4-part query construction
  qp_leg <- queries_progress(
    list(A = "Q123"),
    start = "1990",
    end = "2010",
    limit = "25",
    query_part_1 = "P1:",
    query_part_2 = ":P2:",
    query_part_3 = ":P3:",
    query_part_4 = ":P4"
  )
  expect_identical(qp_leg$A, "P1:Q123:P2:1990:P3:2010:P4\nLIMIT 25")
})

# --- 6. tables_progress (non-empty and empty branches) ---
sc_mock <- data.table::data.table(
  structure_id = c("IK1", "IK2"),
  chemical_pathway = c("Terpenoids", "Alkaloids"),
  chemical_superclass = c("Monoterpenoids", "Indole alkaloids"),
  chemical_class = c("Secoiridoids", "Simple indole")
)
wiki_nonempty <- data.table::data.table(
  structure = c("Q1", "Q2"),
  structureLabel = c("CompA", "CompB"),
  structure_id = c("IK1", "IK2"),
  structureSmiles = c("C=C", NA_character_),
  taxaLabels = c("Swertia chirayita", "Swertia perennis"),
  taxa = c("Q10", "Q11"),
  referencesLabels = c("Ref1", "Ref2"),
  references_ids = c("10.1/1", "10.1/2"),
  references = c("Q100", "Q101")
)
wiki_empty <- data.table::data.table()

tp_res <- tables_progress(
  list(nonempty = wiki_nonempty, empty = wiki_empty),
  structures_classified = sc_mock
)
expect_identical(nrow(tp_res$nonempty), 2L)
expect_false(any(c("references_ids", "structure_id", "structureSmiles") %in% names(tp_res$nonempty)))
expect_identical( names(tp_res$nonempty)[1:3], c("structureLabel", "structure", "structureImage"))
expect_true(any(is.na(tp_res$nonempty$structureImage)))
expect_true(any(tp_res$nonempty$structureImage == utils::URLencode("C=C")))

# Empty table branch returns 1-row NA template
expect_identical(nrow(tp_res$empty), 1L)
expect_true(all(is.na(tp_res$empty$structure)))
expect_identical(
  names(tp_res$empty),
  c(
    "structureLabel",
    "structure",
    "structureImage",
    "taxaLabels",
    "taxa",
    "referencesLabels",
    "references",
    "chemical_pathway",
    "chemical_superclass",
    "chemical_class"
  )
)

# --- 7. format_gt branches (NA/empty link values, http:// URLs vs QIDs, missing structureImage) ---
df_gt_edge <- data.frame(
  structureLabel = c("CompA", "CompB"),
  structure = c("https://www.wikidata.org/entity/Q1", NA_character_),
  taxaLabels = c("Swertia", "Gentiana"),
  taxa = c("Q158572", ""),
  referencesLabels = c("Ref1", "Ref2"),
  references = c("Q999", "http://example.org/ref"),
  stringsAsFactors = FALSE
)
gt_no_img <- format_gt(df_gt_edge, title = "Title", subtitle = "Sub")
expect_inherits(gt_no_img, "gt_tbl")
