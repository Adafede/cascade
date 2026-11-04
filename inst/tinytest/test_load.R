library(tinytest)
library(cascade)

check_export_dir <- utils::getFromNamespace("check_export_dir", "cascade")
fread_gz <- utils::getFromNamespace("fread_gz", "cascade")
get_zenodo_file <- utils::getFromNamespace("get_zenodo_file", "cascade")
load_annotations <- utils::getFromNamespace("load_annotations", "cascade")
load_chromatograms <- utils::getFromNamespace("load_chromatograms", "cascade")
load_features <- utils::getFromNamespace("load_features", "cascade")
load_features_informed <- utils::getFromNamespace(
  "load_features_informed",
  "cascade"
)
load_features_not_informed <- utils::getFromNamespace(
  "load_features_not_informed",
  "cascade"
)
load_ms_data <- utils::getFromNamespace("load_ms_data", "cascade")
parse_ms_spectra <- utils::getFromNamespace("parse_ms_spectra", "cascade")
load_name <- utils::getFromNamespace("load_name", "cascade")
cascade_env <- utils::getFromNamespace(".cascade_env", "cascade")

work_dir <- tempfile("cascade_load_test_")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

withr::with_dir(work_dir, {
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  # --- 1. check_export_dir ---
  nested_dir <- file.path(work_dir, "a", "b", "c")
  expect_false(dir.exists(nested_dir))
  ret_dir <- check_export_dir(nested_dir)
  expect_identical(ret_dir, nested_dir)
  expect_true(dir.exists(nested_dir))
  # Idempotent when directory already exists
  expect_identical(check_export_dir(nested_dir), nested_dir)

  # --- 2. fread_gz (plain file, .gz cold cache, .gz warm cache, .gz evicted
  # file) ---
  plain_csv <- file.path(work_dir, "plain.csv")
  writeLines(c("id,val", "1,10.5", "2,20.5"), plain_csv)
  dt_plain <- fread_gz(plain_csv)
  expect_identical(dt_plain$id, c(1L, 2L))
  expect_equal(dt_plain$val, c(10.5, 20.5), tolerance = 1e-12)

  gz_path <- file.path(work_dir, "sample.csv.gz")
  gz_con <- gzfile(gz_path, open = "wb")
  writeLines(c("id,val", "3,30.25", "4,40.75"), gz_con)
  close(gz_con)

  dt_gz1 <- fread_gz(gz_path)
  expect_identical(dt_gz1$id, c(3L, 4L))
  expect_equal(dt_gz1$val, c(30.25, 40.75), tolerance = 1e-12)
  # Second read hits session cache
  dt_gz2 <- fread_gz(gz_path)
  expect_identical(dt_gz2, dt_gz1)
  # Evict cached decompressed file and verify transparent re-extraction
  cache_key <- paste0("gz_", normalizePath(gz_path, mustWork = FALSE))
  unlink(cascade_env[[cache_key]])
  dt_gz3 <- fread_gz(gz_path)
  expect_identical(dt_gz3, dt_gz1)

  # --- 3. get_zenodo_file branches ---
  orig_session_lotus <- cascade_env[["zenodo_lotus_path"]]
  on.exit(
    {
      cascade_env[["zenodo_lotus_path"]] <- orig_session_lotus
    },
    add = TRUE
  )

  # 3a. Existing non-empty file at destination returns immediately
  existing_dest <- file.path(work_dir, "existing_lotus.csv.gz")
  writeLines("dummy_content", existing_dest)
  expect_identical(get_zenodo_file(path = existing_dest), existing_dest)

  # 3b. Session cache hit (.cascade_env[["zenodo_lotus_path"]])
  cascade_env[["zenodo_lotus_path"]] <- existing_dest
  session_dest <- file.path(work_dir, "sub_session", "from_session.csv.gz")
  expect_identical(get_zenodo_file(path = session_dest), session_dest)
  expect_true(file.exists(session_dest))
  expect_identical(readLines(session_dest), "dummy_content")

  # 3c. User cache directory hit (tools::R_user_dir("cascade", which = "cache"))
  cascade_env[["zenodo_lotus_path"]] <- NULL
  user_cache_dir <- tools::R_user_dir("cascade", which = "cache")
  dir.create(user_cache_dir, recursive = TRUE, showWarnings = FALSE)
  unique_name <- paste0("test_user_cache_", Sys.getpid(), ".csv.gz")
  user_cache_file <- file.path(user_cache_dir, unique_name)
  writeLines("user_cached_content", user_cache_file)
  on.exit(unlink(user_cache_file), add = TRUE)

  uc_dest <- file.path(work_dir, "sub_uc", unique_name)
  expect_identical(get_zenodo_file(path = uc_dest), uc_dest)
  expect_true(file.exists(uc_dest))
  expect_identical(readLines(uc_dest), "user_cached_content")

  # 3d. Zenodo API download & error paths (isolated via mocked curl bindings)
  cascade_env[["zenodo_lotus_path"]] <- NULL
  curl_ns <- asNamespace("curl")
  orig_fetch <- get("curl_fetch_memory", envir = curl_ns)
  orig_dl <- get("curl_download", envir = curl_ns)
  unlockBinding("curl_fetch_memory", curl_ns)
  unlockBinding("curl_download", curl_ns)
  on.exit(
    {
      assign("curl_fetch_memory", orig_fetch, envir = curl_ns)
      assign("curl_download", orig_dl, envir = curl_ns)
      lockBinding("curl_fetch_memory", curl_ns)
      lockBinding("curl_download", curl_ns)
    },
    add = TRUE
  )

  # HTTP 500 error branch
  assign(
    "curl_fetch_memory",
    function(url, handle = NULL) {
      list(status_code = 500L, content = charToRaw("Internal Server Error"))
    },
    envir = curl_ns
  )
  expect_error(
    get_zenodo_file(path = file.path(work_dir, "err500.csv.gz")),
    pattern = "Zenodo API request failed \\(HTTP 500\\)"
  )

  # Pattern not found in Zenodo record JSON branch
  assign(
    "curl_fetch_memory",
    function(url, handle = NULL) {
      list(
        status_code = 200L,
        content = charToRaw(paste0(
          '{"files": [{"links": {"self": ',
          '"https://zenodo.org/api/records/123/files/other.txt"}}]}'
        ))
      )
    },
    envir = curl_ns
  )
  expect_error(
    get_zenodo_file(
      pattern = "frozen_metadata.csv.gz",
      path = file.path(work_dir, "notfound.csv.gz")
    ),
    pattern = "Could not find file matching 'frozen_metadata.csv.gz'"
  )

  # Successful Zenodo JSON match + curl_download
  assign(
    "curl_fetch_memory",
    function(url, handle = NULL) {
      list(
        status_code = 200L,
        content = charToRaw(paste0(
          '{"links": "https://zenodo.org/api/records/5794106/',
          'files/230101_frozen_metadata.csv.gz/content"}'
        ))
      )
    },
    envir = curl_ns
  )
  assign(
    "curl_download",
    function(url, destfile, quiet = TRUE, mode = "wb", handle = NULL) {
      writeLines(paste("downloaded_from", url), destfile)
      invisible(destfile)
    },
    envir = curl_ns
  )
  dl_dest <- file.path(work_dir, "dl", "zenodo_dl.csv.gz")
  expect_identical(get_zenodo_file(path = dl_dest), dl_dest)
  expect_identical(
    readLines(dl_dest),
    paste0(
      "downloaded_from https://zenodo.org/api/records/5794106/",
      "files/230101_frozen_metadata.csv.gz/content"
    )
  )

  # Restore curl bindings immediately
  assign("curl_fetch_memory", orig_fetch, envir = curl_ns)
  assign("curl_download", orig_dl, envir = curl_ns)
  cascade_env[["zenodo_lotus_path"]] <- orig_session_lotus

  # --- 4. load_annotations ---
  ann_file <- file.path(work_dir, "custom_annotations.tsv")
  writeLines(
    c(
      "feature_id\tcandidate_structure_tax_npc_01pat\tother_col",
      "1\tAlkaloids$Terpenoidsu00A7Sub\tkeep$u00A7"
    ),
    ann_file
  )
  dt_ann <- load_annotations(ann_file, show_example = FALSE, mode = "neg")
  expect_identical(
    dt_ann$candidate_structure_tax_npc_01pat,
    "AlkaloidsorTerpenoids$Sub"
  )
  expect_identical(dt_ann$other_col, "keep$u00A7")
  expect_identical(dt_ann$mode, "neg")

  # --- 5. load_name ---
  expect_identical(
    load_name("/some/dir/sample_run.mzML", show_example = FALSE),
    "sample_run.mzML"
  )
  expect_identical(
    load_name(
      "/some/dir/sample_run.mzML",
      default = "def.mzML",
      show_example = TRUE
    ),
    "def.mzML"
  )

  # --- 6. load_features, load_features_informed, load_features_not_informed ---
  feat_tsv <- file.path(work_dir, "feat.tsv")
  writeLines(c("id\trt\tmz", "10\t5.5\t250.1"), feat_tsv)
  dt_feat <- load_features(feat_tsv, show_example = FALSE)
  expect_identical(dt_feat$id, 10L)

  # Fallback to package data when default interim TSV does not exist
  dt_inf_pkg <- load_features_informed(show_example = TRUE)
  expect_identical(nrow(dt_inf_pkg), 3860L)
  dt_ninf_pkg <- load_features_not_informed(show_example = TRUE)
  expect_identical(nrow(dt_ninf_pkg), 893L)

  # Explicit file path when show_example = FALSE and show_example = TRUE
  inf_custom <- file.path(work_dir, "inf_custom.tsv")
  writeLines(c("peak_id\tfeature_id", "1\t99"), inf_custom)
  expect_identical(
    load_features_informed(inf_custom, show_example = FALSE)$feature_id,
    99L
  )
  expect_identical(
    load_features_informed(inf_custom, show_example = TRUE)$feature_id,
    99L
  )
  expect_identical(
    load_features_not_informed(inf_custom, show_example = FALSE)$feature_id,
    99L
  )
  expect_identical(
    load_features_not_informed(inf_custom, show_example = TRUE)$feature_id,
    99L
  )

  # Fallback when file is non-existent and show_example = FALSE (calls
  # fread_gz(file))
  expect_identical(
    load_features_informed(file = inf_custom, show_example = FALSE)$peak_id,
    1L
  )

  # --- 7. parse_ms_spectra & load_ms_data & load_chromatograms ---
  # 7a. Legacy S4/environment assayData with unsorted scans
  legacy_env <- new.env(parent = emptyenv())
  s1 <- list()
  attr(s1, "rt") <- 30.0
  attr(s1, "mz") <- c(100.1, 200.2)
  attr(s1, "intensity") <- c(500, 1000)
  s2 <- list()
  attr(s2, "rt") <- 10.0
  attr(s2, "mz") <- c(150.1)
  attr(s2, "intensity") <- c(250)
  legacy_env[["scan_b"]] <- s1
  legacy_env[["scan_a"]] <- s2
  legacy_obj <- structure(list(), assayData = legacy_env)

  parsed_leg <- parse_ms_spectra(legacy_obj)
  expect_identical(parsed_leg$rt, c(10.0, 30.0))
  expect_identical(parsed_leg$mz, list(c(150.1), c(100.1, 200.2)))
  expect_identical(parsed_leg$intensity, list(c(250), c(500, 1000)))

  # Invalid ms_data object throws informative error
  expect_error(
    parse_ms_spectra(list(wrong = 1)),
    pattern = "Unsupported `ms_data` object format"
  )

  # 7b. Cold cache vs warm cache in load_ms_data(show_example = TRUE)
  orig_ms_cache <- cascade_env[["example_ms_spectra"]]
  on.exit(
    {
      cascade_env[["example_ms_spectra"]] <- orig_ms_cache
    },
    add = TRUE
  )
  cascade_env[["example_ms_spectra"]] <- NULL
  ms_cold <- load_ms_data(show_example = TRUE)
  expect_true(is.list(ms_cold))
  expect_identical(names(ms_cold), c("rt", "mz", "intensity"))
  ms_warm <- load_ms_data(show_example = TRUE)
  expect_identical(ms_warm$rt, ms_cold$rt)

  # 7c. Real mzML file loading for load_ms_data and load_chromatograms
  # Base64-encoded 64-bit little-endian IEEE754 floats for 100/200, 10/20, 1/2:
  mzml_1scan <- paste0(
    "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
    "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\n",
    "  <run id=\"r1\">\n",
    "    <spectrumList count=\"1\">\n",
    "      <spectrum index=\"0\" id=\"scan=1\" defaultArrayLength=\"2\">\n",
    "        <cvParam cvRef=\"MS\" accession=\"MS:1000511\"",
    " name=\"ms level\" value=\"1\"/>\n",
    "        <scanList count=\"1\"><scan>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000016\"",
    " name=\"scan start time\" value=\"15.0\"",
    " unitCvRef=\"UO\" unitAccession=\"UO:0000010\"",
    " unitName=\"second\"/></scan></scanList>\n",
    "        <binaryDataArrayList count=\"2\">\n",
    "          <binaryDataArray encodedLength=\"24\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"/>",
    "<binary>AAAAAAAAWUAAAAAAAABpQA==</binary></binaryDataArray>\n",
    "          <binaryDataArray encodedLength=\"24\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>",
    "<binary>AAAAAAAAJEAAAAAAAAA0QA==</binary></binaryDataArray>\n",
    "        </binaryDataArrayList>\n",
    "      </spectrum>\n",
    "    </spectrumList>\n",
    "    <chromatogramList count=\"1\">\n",
    "      <chromatogram index=\"0\" id=\"UV#1_CAD_1_0\"",
    " defaultArrayLength=\"2\">\n",
    "        <binaryDataArrayList count=\"2\">\n",
    "          <binaryDataArray encodedLength=\"24\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000595\" name=\"time array\"/>",
    "<binary>AAAAAAAA8D8AAAAAAAAAQA==</binary></binaryDataArray>\n",
    "          <binaryDataArray encodedLength=\"24\">",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\"/>",
    "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>",
    "<binary>AAAAAAAAJEAAAAAAAAA0QA==</binary></binaryDataArray>\n",
    "        </binaryDataArrayList>\n",
    "      </chromatogram>\n",
    "    </chromatogramList>\n",
    "  </run>\n",
    "</mzML>"
  )
  mzml_file <- file.path(work_dir, "single.mzML")
  writeLines(mzml_1scan, mzml_file)

  ms_from_file <- load_ms_data(mzml_file, show_example = FALSE)
  expect_equal(ms_from_file$rt, 15.0, tolerance = 1e-12)
  expect_equal(ms_from_file$mz[[1L]], c(100, 200), tolerance = 1e-12)
  expect_equal(ms_from_file$intensity[[1L]], c(10, 20), tolerance = 1e-12)

  chrom_from_file <- load_chromatograms(mzml_file, show_example = FALSE)
  expect_identical(names(chrom_from_file), "cad")
  expect_equal(chrom_from_file$cad$rtime, c(1, 2), tolerance = 1e-12)
  expect_equal(chrom_from_file$cad$intensity, c(10, 20), tolerance = 1e-12)

  # Negative polarity example chromatograms
  chrom_neg <- load_chromatograms(show_example = TRUE, example_polarity = "neg")
  expect_identical(names(chrom_neg), c("bpi", "pda", "cad"))
  expect_true(nrow(chrom_neg$cad) > 0L)
})
