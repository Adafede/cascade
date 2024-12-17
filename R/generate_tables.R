#' Generate IDs
#'
#' @export
#'
#' @include check_export_dir.R
#' @include keep_best_candidates.R
#' @include load_annotations.R
#' @include load_features_informed.R
#' @include make_confident.R
#' @include molinfo.R
#'
#' @param annotations Annotations
#' @param file_negative File negative
#' @param file_positive File positive
#' @param min_confidence Min confidence
#' @param show_example Show example? Default to FALSE
#' @param export_csv Export CSV
#' @param export_html Export HTML
#' @param export_dir Export Dir
#' @param export_name Export name
#'
#' @return Tables
#'
#' @examples
#' \dontrun{
#' generate_tables()
#' }
generate_tables <- function(annotations = NULL,
                            file_negative = NULL,
                            file_positive = NULL,
                            min_confidence = 0.4,
                            show_example = FALSE,
                            export_csv = TRUE,
                            export_html = TRUE,
                            export_dir = "data/processed",
                            export_name = "cascade_table") {
  message("loading informed features")
  if (!is.null(file_positive) | show_example) {
    message("... positive mode")
    if (show_example) {
      tables_peaks_pos <- load_features_informed(show_example = show_example)
    } else {
      tables_peaks_pos <- file_positive |>
        load_features_informed()
    }
  }
  if (!is.null(file_negative)) {
    message("... negative mode")
    tables_peaks_neg <- file_negative |>
      load_features_informed(show_example = show_example) |>
      tidytable::select(tidytable::all_of(
        c("peak_rt_min", "peak_rt_apex", "peak_rt_max", "peak_area")
      )) |>
      tidytable::distinct()
  }

  message("loading annotations")
  annotation_table <- annotations |>
    load_annotations(show_example = show_example) |>
    keep_best_candidates() |>
    tidytable::mutate(feature_id = as.numeric(feature_id)) |>
    make_confident(score = min_confidence)

  message("Getting last LOTUS version")
  tima::get_last_version_from_zenodo(
    doi = "10.5281/zenodo.5794106",
    pattern = "frozen_metadata.csv.gz",
    "data/source/libraries/lotus.csv.gz"
  )

  message("Loading LOTUS classified structures")
  lotus <- tidytable::fread(
    file = "data/source/libraries/lotus.csv.gz",
    select = c(
      "structure_inchikey",
      "organism_taxonomy_01domain",
      "organism_taxonomy_02kingdom",
      "organism_taxonomy_03phylum",
      "organism_taxonomy_04class",
      "organism_taxonomy_05order",
      "organism_taxonomy_06family",
      "organism_taxonomy_07tribe",
      "organism_taxonomy_08genus",
      "organism_taxonomy_09species",
      "organism_taxonomy_10varietas",
      "reference_doi"
    )
  ) |>
    tidytable::distinct()


  table_peaks <- tables_peaks_pos

  if (exists("tables_peaks_neg")) {
    table_peaks <- table_peaks |>
      dplyr::bind_rows(table_peaks_neg)
  }

  data.table::setkey(table_peaks, peak_rt_min, peak_rt_max)

  table_processed <- table_peaks |>
    tidytable::inner_join(annotation_table) |>
    tidytable::select(
      peak_rt_min,
      peak_rt_max,
      feature_id,
      feature_rt,
      feature_mz,
      smiles_no_stereo = smiles_2D,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      score_final,
      best_candidate_organism,
      inchikey_no_stereo = inchikey_2D
    )

  data.table::setkey(table_processed, peak_rt_min, peak_rt_max)

  message("joining within given rt tolerance")
  table_medium <-
    data.table::foverlaps(table_peaks, table_processed) |>
    tidytable::arrange(tidytable::desc(score_final)) |>
    tidytable::distinct(peak_rt_apex, smiles_no_stereo, .keep_all = TRUE) |>
    tidytable::group_by(peak_rt_min) |>
    tidytable::mutate(
      newrt = mean(peak_rt_apex),
      newarea = mean(peak_area)
    ) |>
    tidytable::ungroup() |>
    tidytable::distinct(newrt, newarea, .keep_all = TRUE) |>
    tidytable::arrange(newrt) |>
    tidytable::mutate(diff = newrt - tidytable::lag(newrt)) |> # eliminate rounding duplicates
    tidytable::filter(!is.na(diff) &
      diff > 0.05) |>
    tidytable::mutate(peak_id = tidytable::row_number()) |>
    tidytable::select(
      peak_id,
      peak_rt = newrt,
      peak_area = newarea,
      feature_id,
      feature_rt,
      feature_mz,
      smiles_no_stereo,
      best_candidate_1,
      best_candidate_2,
      best_candidate_3,
      score_final,
      best_candidate_organism,
      inchikey_no_stereo
    ) |>
    tidytable::mutate(
      peak_rt = round(x = peak_rt, digits = 2),
      peak_area = round(x = 100 * peak_area / sum(peak_area), digits = 2),
      feature_rt = round(x = feature_rt, digits = 2),
      feature_mz = round(x = feature_mz, digits = 4),
      score_final = round(x = as.numeric(score_final), digits = 2)
    )

  lotus_long <- lotus |>
    tidytable::mutate(inchikey_no_stereo = gsub(
      pattern = "-.*",
      replacement = "",
      x = structure_inchikey
    )) |>
    tidytable::pivot_longer(cols = 2:11, names_prefix = "organism_taxonomy_") |>
    tidytable::select(inchikey_no_stereo,
      reference = reference_doi,
      best_candidate_organism = value
    ) |>
    tidytable::distinct()

  table_final <- table_medium |>
    tidytable::left_join(lotus_long) |>
    tidytable::group_by(peak_id, inchikey_no_stereo) |>
    tidytable::mutate(reference = paste(reference, collapse = "; ")) |>
    tidytable::mutate(Structure = URLencode(smiles_no_stereo)) |>
    tidytable::mutate(
      reference = ifelse(
        test = reference == "NA",
        yes = NA,
        no = reference
      ),
      Structure = ifelse(
        test = Structure == "NA",
        yes = "",
        no = Structure
      ),
      score_final = ifelse(
        test = score_final > 1,
        yes = as.character(">1.00"),
        no = as.character(score_final)
      )
    ) |>
    tidytable::ungroup() |>
    tidytable::distinct() |>
    tidytable::select(
      `Peak ID` = peak_id,
      `Peak RT [min]` = peak_rt,
      `Peak Area [%]` = peak_area,
      `Feature ID` = feature_id,
      `Feature RT [min]` = feature_rt,
      `Feature m/z` = feature_mz,
      Structure,
      `InChIKey no stereo` = inchikey_no_stereo,
      `Annotation Score` = score_final,
      `Chemical Pathway` = best_candidate_1,
      `Chemical Superclass` = best_candidate_2,
      `Chemical Class` = best_candidate_3,
      `Closest Organism with Structure Reported` = best_candidate_organism,
      `Reference(s)` = reference
    )

  table_pretty <- table_final |>
    gt::gt(auto_align = TRUE) |>
    # cols_width(
    #   everything() ~ px(200)
    # ) |>
    gt::text_transform(
      locations = gt::cells_body(columns = Structure),
      fn = molinfo
    ) |>
    gt::sub_missing(
      columns = tidytable::everything(),
      rows = tidytable::everything(),
      missing_text = ""
    ) |>
    gt::opt_interactive(use_filters = TRUE)

  if (export_csv | export_html) {
    message("checking export directory")
    check_export_dir(export_dir)
    message("exporting")
  }

  if (export_csv) {
    message("...CSV")
    table_final |>
      tidytable::fwrite(file = file.path(export_dir, paste0(export_name, ".tsv")), sep = "\t")
  }
  if (export_html) {
    message("...HTML")
    table_pretty |>
      gt::gtsave(filename = file.path(export_dir, paste0(export_name, ".html")))
  }

  return(list("pretty_table" = table_pretty, "csv_table" = table_final))
}
