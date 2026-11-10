#' CASCADE Example Datasets
#'
#' Example datasets included with CASCADE for demonstration and testing
#' purposes. All tabular datasets are stored as [data.frame] objects and
#' converted to [data.table::data.table] when loaded via the internal `load_*`
#' helpers.
#'
#' @section cascade_annotations:
#' Metabolite annotation table from a *Swertia* extract analyzed with TIMA.
#' Contains 3,471 unique features with their chemical and taxonomic
#' contextualizations.
#'
#' @section cascade_features:
#' Complete mzmine feature table. Contains 3,471 features with *m/z*,
#' retention time, and intensity.
#'
#' @section cascade_features_informed:
#' CAD-detected peaks that matched annotation hypotheses. Contains 4,132
#' peak-feature matches.
#'
#' @section cascade_features_not_informed:
#' CAD-detected peaks that did not match any annotation hypothesis. Contains
#' 650 peak-feature matches.
#'
#' @section cascade_chromatograms_positive:
#' Chromatogram traces from positive-mode analysis. Tidy long-form
#' [data.frame] with 3 columns.
#'
#' @section cascade_chromatograms_negative:
#' Chromatogram traces from negative-mode analysis. Tidy long-form
#' [data.frame] with 3 columns.
#'
#' @section cascade_ms_data:
#' Example MS1 spectral dataset used for extracted ion chromatogram (XIC)
#' reconstruction and peak shape comparison.
#'
#' @format
#' - **cascade_annotations**: [data.frame], 3,471 rows x 58 columns
#'   (EDAM: `data_0962` Small molecule report)
#' - **cascade_features**: [data.frame], 3,471 rows x 101 columns
#'   (EDAM: `data_2536` Mass spectrometry data)
#' - **cascade_features_informed**: [data.frame], 4,132 rows x 11 columns
#'   (EDAM: `data_2536` Mass spectrometry data)
#' - **cascade_features_not_informed**: [data.frame], 650 rows x 11 columns
#'   (EDAM: `data_2536` Mass spectrometry data)
#' - **cascade_chromatograms_positive**: [data.frame], 122,541 rows x 3 columns
#'   (EDAM: `data_2536` Mass spectrometry data)
#' - **cascade_chromatograms_negative**: [data.frame], 122,079 rows x 3 columns
#'   (EDAM: `data_2536` Mass spectrometry data)
#' - **cascade_ms_data**: Named list (`rt`, `mz`, `intensity`) containing
#'   5,668 MS1 scans (EDAM: `data_0943` Mass spectrum)
#'
#' @concept EDAM:data_2536 Mass spectrometry data
#' @concept EDAM:data_0943 Mass spectrum
#' @concept EDAM:data_0962 Small molecule report
#' @concept EDAM:topic_3172 Metabolomics
#' @concept EDAM:topic_0154 Small molecules
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
