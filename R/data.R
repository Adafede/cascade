#' CASCADE Example Datasets
#'
#' Example datasets included with CASCADE for demonstration and testing purposes.
#' All data is stored as [data.frame] format for efficiency, and is converted to
#' [tidytable::tidytable] when loaded via the load_* functions.
#'
#' @section cascade_annotations:
#' Metabolite annotation table from a Swertia extract analyzed with TIMA.
#' Contains 3471 unique features with their chemical and taxonomic contextualizations.
#'
#' @section cascade_features:
#' Complete mzMine feature table. Contains 3471 features with m/z, retention time, and intensity.
#'
#' @section cascade_features_informed:
#' CAD-detected peaks that matched annotation hypotheses. Contains 3800 peak-feature matches.
#'
#' @section cascade_features_not_informed:
#' CAD-detected peaks that did NOT match any annotation hypothesis. Contains 910 peak-feature matches.
#'
#' @section cascade_chromatograms_positive:
#' Chromatogram traces from positive-mode analysis. Tidy long-form data.frame with 3 columns.
#'
#' @section cascade_chromatograms_negative:
#' Chromatogram traces from negative-mode analysis. Tidy long-form data.frame with 3 columns.
#'
#' @section cascade_ms_data:
#' Example MS data as an [MSnbase::MSnExp] object. Required for full MS workflow.
#'
#' @format
#'
#' **cascade_annotations**: [data.frame], 3471 rows × 58 columns
#' **cascade_features**: [data.frame], 3471 rows × 101 columns
#' **cascade_features_informed**: [data.frame], 3800 rows × 11 columns
#' **cascade_features_not_informed**: [data.frame], 910 rows × 11 columns
#' **cascade_chromatograms_positive**: [data.frame], 122541 rows × 3 columns
#' **cascade_chromatograms_negative**: [data.frame], 122079 rows × 3 columns
#' **cascade_ms_data**: [MSnbase::MSnExp] object
#'
#' @keywords datasets
#' @name cascade_datasets
NULL

#' @rdname cascade_datasets
"cascade_annotations"

#' @rdname cascade_datasets
"cascade_features"

#' @rdname cascade_datasets
"cascade_features_informed"

#' @rdname cascade_datasets
"cascade_features_not_informed"

#' @rdname cascade_datasets
"cascade_chromatograms_positive"

#' @rdname cascade_datasets
"cascade_chromatograms_negative"

#' @rdname cascade_datasets
"cascade_ms_data"
