#' CASCADE Example Datasets
#'
#' This document describes the example datasets included with CASCADE for
#' demonstration and testing purposes.
#'
#' @section cascade_annotations:
#' Metabolite annotation table from a Swertia extract analyzed with TIMA
#' (Taxonomically Informed Microbial Annotation).
#' Contains 3471 unique features with their chemical and taxonomic
#' contextualizations.
#'
#' @section cascade_features:
#' Complete mzMine feature table for the example extract.
#' Contains 3471 features with m/z, retention time, and intensity information.
#'
#' @section cascade_features_informed:
#' Filtered feature table containing only CAD-detected peaks that matched
#' annotation hypotheses from the TIMA table.
#' Contains 3800 peak-feature matches with semi-quantitative CAD signal integration.
#'
#' @section cascade_features_not_informed:
#' Filtered feature table containing only CAD-detected peaks that did NOT match
#' any annotation hypothesis. Contains 910 peak-feature matches.
#'
#' @section cascade_chromatograms_positive:
#' Tidy long-form chromatogram traces from the positive-mode example mzML file.
#' Stored as a data.frame with columns `chromatogram`, `rtime`, and `intensity`.
#'
#' @section cascade_chromatograms_negative:
#' Tidy long-form chromatogram traces from the negative-mode example mzML file.
#' Stored as a data.frame with columns `chromatogram`, `rtime`, and `intensity`.
#'
#' @section cascade_ms_data:
#' Example MS data used for chromatogram extraction and peak workflows.
#' Stored as an [MSnbase::MSnExp] object because the downstream workflow
#' requires the full MS experiment structure.
#'
#' @format
#'
#' **cascade_annotations**: A [data.frame] with 3471 rows and 58 columns
#' including feature identifiers, chemical structure information (InChI, SMILES),
#' and taxonomic contextualizations (NCBI taxonomy ranks).
#'
#' **cascade_features**: A [data.frame] with 3471 rows and 101 columns
#' including m/z, retention time, intensity, peak shape parameters, and manual
#' annotation fields from mzMine.
#'
#' **cascade_features_informed**: A [data.frame] with 3800 rows and 11 columns
#' including sample identifier, peak ID, retention time boundaries (rt_min, rt_apex, rt_max),
#' peak apex intensity, and CAD-integrated signal (integral column).
#'
#' **cascade_features_not_informed**: A [data.frame] with 910 rows and 11 columns
#' with the same schema as cascade_features_informed but for unannotated peaks.
#'
#' **cascade_chromatograms_positive**: A [data.frame] with 122541 rows and 3 columns
#' (`chromatogram`, `rtime`, `intensity`) containing positive-mode example traces.
#'
#' **cascade_chromatograms_negative**: A [data.frame] with 122079 rows and 3 columns
#' (`chromatogram`, `rtime`, `intensity`) containing negative-mode example traces.
#'
#' **cascade_ms_data**: An [MSnbase::MSnExp] object containing the example MS data
#' used by the workflow.
#'
#' @source
#' All datasets are derived from a Swertia chirayita (Roxb.) H.Karst. extract
#' analyzed by DDA-MS with CAD detection as described in:
#'
#' Rutz & Wolfender (2023) "Automated Composition Assessment of Natural Extracts:
#' Untargeted Mass Spectrometry-Based Metabolite Profiling Integrating
#' Semiquantitative Detection" \emph{Journal of Agricultural and Food Chemistry}
#' 71(46):18010-18023. \doi{10.1021/acs.jafc.3c03099}
#'
#' @usage
#' cascade_annotations
#' cascade_features
#' cascade_features_informed
#' cascade_features_not_informed
#' cascade_chromatograms_positive
#' cascade_chromatograms_negative
#' cascade_ms_data
#'
#' @keywords datasets
#'
#' @examples
#' # Load CASCADE example datasets
#' data(cascade_annotations)
#' data(cascade_features)
#' data(cascade_features_informed)
#' data(cascade_features_not_informed)
#' data(cascade_chromatograms_positive)
#' data(cascade_chromatograms_negative)
#' data(cascade_ms_data)
#'
#' # Check dimensions
#' dim(cascade_annotations)
#' dim(cascade_features_informed)
#'
#' # Explore column names
#' names(cascade_features)
#'
#' @name cascade_datasets
#' @aliases cascade_annotations cascade_features cascade_features_informed cascade_features_not_informed cascade_chromatograms_positive cascade_chromatograms_negative cascade_ms_data
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
