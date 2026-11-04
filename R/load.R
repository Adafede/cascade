#' Ensure an export directory exists
#'
#' @param dir Character(1) directory path to create if missing.
#'
#' @return Invisibly returns `dir`.
#' @keywords internal
#' @noRd
check_export_dir <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dir)
}

#' Read a delimited file, transparently decompressing .gz via base R zlib
#'
#' @param file Character(1) path to a delimited text or `.gz` file.
#' @param ... Additional arguments passed to [data.table::fread()].
#'
#' @return A `data.table`.
#' @keywords internal
#' @noRd
fread_gz <- function(file, ...) {
  if (grepl("\\.gz$", file, ignore.case = TRUE)) {
    cache_key <- paste0("gz_", normalizePath(file, mustWork = FALSE))
    cached_tmp <- .cascade_env[[cache_key]]
    if (is.null(cached_tmp) || !file.exists(cached_tmp)) {
      cached_tmp <- tempfile(pattern = "cascade_gz_", fileext = ".csv")
      con_in <- gzfile(file, open = "rb")
      on.exit(close(con_in), add = TRUE)
      con_out <- file(cached_tmp, open = "wb")
      on.exit(close(con_out), add = TRUE)
      repeat {
        buf <- readBin(con_in, what = "raw", n = 16777216L)
        if (length(buf) == 0L) {
          break
        }
        writeBin(buf, con_out)
      }
      close(con_out)
      close(con_in)
      on.exit(NULL, add = FALSE)
      .cascade_env[[cache_key]] <- cached_tmp
    }
    return(data.table::fread(file = cached_tmp, showProgress = FALSE, ...))
  }
  data.table::fread(file = file, showProgress = FALSE, ...)
}

#' Download or locate the latest matching file from a Zenodo record DOI
#'
#' @param doi Character(1) Zenodo DOI (default `"10.5281/zenodo.5794106"`).
#' @param pattern Character(1) filename suffix pattern to match in the record.
#' @param path Character(1) destination file path.
#' @param timeout_s Integer(1) download timeout in seconds (default `300L`).
#'
#' @return Invisibly returns `path`.
#' @keywords internal
#' @noRd
get_zenodo_file <- function(
  doi = "10.5281/zenodo.5794106",
  pattern = "frozen_metadata.csv.gz",
  path = "data/source/libraries/lotus.csv.gz",
  timeout_s = 300L
) {
  if (file.exists(path) && file.size(path) > 0) {
    return(invisible(path))
  }

  # Check session cache so repeated calls across tempdirs reuse the file
  session_cached <- .cascade_env[["zenodo_lotus_path"]]
  if (
    !is.null(session_cached) &&
      file.exists(session_cached) &&
      file.size(session_cached) > 0
  ) {
    check_export_dir(dirname(path))
    if (!file.symlink(session_cached, path)) {
      file.copy(session_cached, path, overwrite = TRUE)
    }
    return(invisible(path))
  }

  # Also check user cache directory
  cache_dir <- tools::R_user_dir("cascade", which = "cache")
  cache_file <- file.path(cache_dir, basename(path))
  if (file.exists(cache_file) && file.size(cache_file) > 0) {
    .cascade_env[["zenodo_lotus_path"]] <- normalizePath(cache_file)
    check_export_dir(dirname(path))
    if (!file.symlink(cache_file, path)) {
      file.copy(cache_file, path, overwrite = TRUE)
    }
    return(invisible(path))
  }

  record_id <- sub("^.*zenodo\\.", "", doi)
  api_url <- paste0("https://zenodo.org/api/records/", record_id)

  h <- curl::new_handle()
  curl::handle_setopt(h, timeout = as.integer(timeout_s), followlocation = 1L)
  resp <- curl::curl_fetch_memory(api_url, handle = h)
  if (resp$status_code >= 400L) {
    stop(
      sprintf(
        "Zenodo API request failed (HTTP %d) for %s",
        resp$status_code,
        api_url
      ),
      call. = FALSE
    )
  }

  txt <- rawToChar(resp$content)
  pat_escaped <- gsub("\\.", "[.]", pattern)
  url_regex <- paste0(
    "https://zenodo[.]org/api/records/[0-9]+/files/[^\"]*",
    pat_escaped,
    "(/content)?"
  )
  m <- regmatches(txt, regexpr(url_regex, txt))
  if (length(m) == 0L || !nzchar(m)) {
    stop(
      sprintf(
        "Could not find file matching '%s' in Zenodo record %s",
        pattern,
        doi
      ),
      call. = FALSE
    )
  }

  check_export_dir(dirname(path))
  dl_handle <- curl::new_handle()
  curl::handle_setopt(
    dl_handle,
    timeout = as.integer(timeout_s),
    followlocation = 1L
  )
  curl::curl_download(
    url = m[1L],
    destfile = path,
    quiet = TRUE,
    handle = dl_handle
  )
  .cascade_env[["zenodo_lotus_path"]] <- normalizePath(path, mustWork = FALSE)
  invisible(path)
}

#' Load annotations
#'
#' @param file Character(1) path to a TIMA annotation TSV/CSV file.
#' @param show_example Logical(1); load bundled `cascade_annotations` dataset.
#' @param mode Character(1) ionization mode (`"pos"` or `"neg"`).
#'
#' @return A `data.table` of structural annotations.
#' @keywords internal
#' @noRd
load_annotations <- function(file = NULL, show_example = FALSE, mode = "pos") {
  if (show_example) {
    utils::data(
      "cascade_annotations",
      package = "cascade",
      envir = environment()
    )
    dt <- data.table::as.data.table(cascade_annotations)
    return(dt)
  }

  dt <- fread_gz(file)
  tax_cols <- grepv(
    "candidate_structure_tax",
    names(dt),
    fixed = TRUE
  )
  for (col in tax_cols) {
    v <- dt[[col]]
    if (is.character(v)) {
      v <- gsub("\\$", "or", v)
      v <- gsub("u00A7", "$", v, fixed = TRUE)
      data.table::set(dt, j = col, value = v)
    }
  }
  data.table::set(dt, j = "mode", value = mode)
  dt
}

#' Load chromatograms
#'
#' @param file Character(1) path to an `.mzML` file.
#' @param headers Named character vector mapping detector keys (`"bpi"`,
#'   `"pda"`, `"cad"`) to chromatogram header IDs.
#' @param show_example Logical(1); load bundled example chromatograms.
#' @param example_polarity Character(1) polarity (`"pos"` or `"neg"`) when
#'   `show_example = TRUE`.
#'
#' @return Named list of chromatogram `data.frame`s (`rtime`, `intensity`).
#' @keywords internal
#' @noRd
load_chromatograms <- function(
  file = NULL,
  headers = c(
    "bpi" = "BasePeak_0",
    "pda" = "PDA#1_TotalAbsorbance_0",
    "cad" = "UV#1_CAD_1_0"
  ),
  show_example = FALSE,
  example_polarity = "pos"
) {
  if (show_example) {
    dataset_name <- switch(
      example_polarity,
      "pos" = "cascade_chromatograms_positive",
      "neg" = "cascade_chromatograms_negative"
    )
    utils::data(list = dataset_name, package = "cascade", envir = environment())
    chromatograms_df <- get(dataset_name, envir = environment())
    chrom_col <- chromatograms_df$chromatogram
    rtime_col <- chromatograms_df$rtime
    int_col <- chromatograms_df$intensity
    chromatograms <- lapply(
      names(headers),
      function(chrom_name) {
        sel <- chrom_col == chrom_name
        data.frame(
          rtime = rtime_col[sel],
          intensity = int_col[sel]
        )
      }
    )
    names(chromatograms) <- names(headers)
    return(chromatograms)
  }

  file_pointer <- mzR::openMSfile(file)
  on.exit(mzR::close(file_pointer), add = TRUE)
  file_headers <- mzR::chromatogramHeader(file_pointer)
  matched <- match(headers, file_headers$chromatogramId)
  valid <- !is.na(matched)
  indices <- file_headers$chromatogramIndex[matched[valid]]
  chromatograms <- mzR::chromatograms(file_pointer, chrom = indices)
  if (is.data.frame(chromatograms)) {
    chromatograms <- list(chromatograms)
  }
  names(chromatograms) <- names(headers)[valid]
  chromatograms
}

#' Load features
#'
#' @param file Character(1) path to a feature CSV/TSV file.
#' @param show_example Logical(1); load bundled `cascade_features`.
#'
#' @return A `data.table` of MS features.
#' @keywords internal
#' @noRd
load_features <- function(file = NULL, show_example = FALSE) {
  if (show_example) {
    utils::data("cascade_features", package = "cascade", envir = environment())
    return(data.table::as.data.table(cascade_features))
  }
  fread_gz(file)
}

#' Load informed features
#'
#' @param file Character(1) path to an informed feature TSV/CSV file.
#' @param show_example Logical(1); load bundled `cascade_features_informed`.
#'
#' @return A `data.table` of informed MS features.
#' @keywords internal
#' @noRd
load_features_informed <- function(file = NULL, show_example = FALSE) {
  if (!is.null(file) && file.exists(file)) {
    return(fread_gz(file))
  }
  if (show_example) {
    default_tsv <- file.path(
      "data/interim/peaks",
      "210619_AR_06_V_03_2_01_featuresInformed_cad.tsv"
    )
    if (file.exists(default_tsv) && file.size(default_tsv) > 0) {
      return(fread_gz(default_tsv))
    }
    utils::data(
      "cascade_features_informed",
      package = "cascade",
      envir = environment()
    )
    return(data.table::as.data.table(cascade_features_informed))
  }
  fread_gz(file)
}

#' Load non-informed features
#'
#' @param file Character(1) path to a non-informed feature TSV/CSV file.
#' @param show_example Logical(1); load bundled `cascade_features_not_informed`.
#'
#' @return A `data.table` of non-informed MS features.
#' @keywords internal
#' @noRd
load_features_not_informed <- function(file = NULL, show_example = FALSE) {
  if (!is.null(file) && file.exists(file)) {
    return(fread_gz(file))
  }
  if (show_example) {
    default_tsv <- file.path(
      "data/interim/peaks",
      "210619_AR_06_V_03_2_01_featuresNotInformed_cad.tsv"
    )
    if (file.exists(default_tsv) && file.size(default_tsv) > 0) {
      return(fread_gz(default_tsv))
    }
    utils::data(
      "cascade_features_not_informed",
      package = "cascade",
      envir = environment()
    )
    return(data.table::as.data.table(cascade_features_not_informed))
  }
  fread_gz(file)
}

#' Load MS1 spectral data
#'
#' @param file Character(1) path to an `.mzML` file.
#' @param show_example Logical(1); load bundled `cascade_ms_data`.
#'
#' @return Named list with elements `rt` (numeric vector), `mz` (list of numeric
#'   vectors), and `intensity` (list of numeric vectors).
#' @keywords internal
#' @noRd
load_ms_data <- function(file = NULL, show_example = FALSE) {
  if (show_example) {
    message(
      "Loading example MS data from package data, ",
      "doing it on disk will be more efficient"
    )
    cached <- .cascade_env[["example_ms_spectra"]]
    if (!is.null(cached)) {
      return(cached)
    }
    utils::data(
      list = "cascade_ms_data",
      package = "cascade",
      envir = environment()
    )
    parsed <- parse_ms_spectra(cascade_ms_data)
    .cascade_env[["example_ms_spectra"]] <- parsed
    return(parsed)
  }

  fh <- mzR::openMSfile(file)
  on.exit(mzR::close(fh), add = TRUE)
  hdr <- mzR::header(fh)
  ms1_idx <- which(hdr$msLevel == 1L)
  rts <- hdr$retentionTime[ms1_idx]
  ord <- order(rts)
  ms1_idx <- ms1_idx[ord]
  rts <- rts[ord]
  pks <- mzR::peaks(fh, scans = ms1_idx)
  if (is.matrix(pks)) {
    pks <- list(pks)
  }
  list(
    rt = rts,
    mz = lapply(pks, function(m) m[, 1L]),
    intensity = lapply(pks, function(m) m[, 2L])
  )
}

#' Convert MSnExp or list spectral container into indexed list
#'
#' @param ms_data Either a pre-parsed `list(rt, mz, intensity)` or a legacy S4
#'   `OnDiskMSnExp`/`MSnExp` object with an `assayData` environment.
#'
#' @return Named list with `rt`, `mz`, and `intensity` sorted by retention time.
#' @keywords internal
#' @noRd
parse_ms_spectra <- function(ms_data) {
  if (is.list(ms_data) && all(c("rt", "mz", "intensity") %in% names(ms_data))) {
    return(ms_data)
  }
  env <- attr(ms_data, "assayData")
  if (is.environment(env)) {
    nms <- ls(env)
    rts <- vapply(
      nms,
      function(k) attr(env[[k]], "rt"),
      numeric(1L),
      USE.NAMES = FALSE
    )
    ord <- order(rts)
    nms <- nms[ord]
    rts <- rts[ord]
    mzs <- lapply(nms, function(k) attr(env[[k]], "mz"))
    ints <- lapply(nms, function(k) attr(env[[k]], "intensity"))
    return(list(rt = rts, mz = mzs, intensity = ints))
  }
  stop("Unsupported `ms_data` object format.", call. = FALSE)
}

#' Load sample name
#'
#' @param file Character(1) file path.
#' @param default Character(1) default sample name when `show_example = TRUE`.
#' @param show_example Logical(1); return `default` instead of `basename(file)`.
#'
#' @return Character(1) sample filename.
#' @keywords internal
#' @noRd
load_name <- function(
  file = NULL,
  default = "210619_AR_06_V_03_2_01.mzML",
  show_example = FALSE
) {
  if (!show_example) {
    basename(file)
  } else {
    default
  }
}
