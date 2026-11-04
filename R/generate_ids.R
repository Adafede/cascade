#' Generate chemotaxonomic literature hierarchies and plots from Wikidata/LOTUS
#'
#' Queries Wikidata for natural products reported in the specified `taxa`
#' within a publication year range (`start` to `end`), enriches them with
#' NPClassifier chemical classifications from LOTUS, optionally deduplicates
#' stereoisomers and filters to chromatographiable mass/XLogP ranges, and
#' builds stacked histograms, interactive treemaps, and sunburst plots.
#'
#' @export
#'
#' @param taxa Character vector of scientific taxon names to query. Default is
#'   `c("Swertia", "Kopsia", "Ginkgo")`.
#' @param comparison Optional character vector of two taxon names in `taxa` for
#'   side-by-side comparative visualization. Default is `NULL`.
#' @param no_stereo Logical; collapse stereoisomers to unique 2D structures per
#'   taxon/reference. Default is `TRUE`.
#' @param filter_ms_conditions Logical; filter structures to chromatographiable
#'   mass (`50..1500` Da) and XLogP (`-1..6`) ranges. Default is `TRUE`.
#' @param start Character(1); minimum publication year. Default is `"0"`.
#' @param end Character(1); maximum publication year. Default is `"9999"`.
#' @param limit Character(1); maximum number of SPARQL rows per taxon. Default
#'   is `"1000000"`.
#' @param show_title Logical; include taxon name and unique structure count
#'   titles on interactive treemaps and sunbursts. Default is `FALSE`.
#' @param format_tables Logical; format the returned `tables` element as
#'   interactive [gt::gt()] tables with clickable Wikidata links and 2D
#'   structure depictions. Default is `FALSE`.
#'
#' @return A named list with four elements:
#'   - `plots`: Named list of `ggplot` chemical class histograms.
#'   - `treemaps`: Named list of interactive `plotly` treemap widgets.
#'   - `sunbursts`: Named list of interactive `plotly` sunburst widgets.
#'   - `tables`: Named list of compound-taxon-reference `data.table` objects
#'     (or `gt_tbl` objects when `format_tables = TRUE`).
#'
#' @concept EDAM:operation_0224 Query and retrieval
#' @concept EDAM:operation_2942 Treemap visualisation
#' @concept EDAM:operation_3803 Natural product identification
#' @concept EDAM:topic_3172 Metabolomics
#' @concept EDAM:topic_0154 Small molecules
#'
#' @examples
#' \donttest{
#' old_wd <- setwd(tempdir())
#' try(generate_ids(taxa = "Swertia"))
#' setwd(old_wd)
#' }
generate_ids <- function(
  taxa = c("Swertia", "Kopsia", "Ginkgo"),
  comparison = NULL,
  no_stereo = TRUE,
  filter_ms_conditions = TRUE,
  start = "0",
  end = "9999",
  limit = "1000000",
  show_title = FALSE,
  format_tables = FALSE
) {
  query_part_1 <- "SELECT DISTINCT ?structure ?structureLabel ?structure_id ?structureSmiles ?taxaLabels ?taxa ?referencesLabels ?references_ids ?references WHERE {\n  ?taxa (wdt:P171*) wd:"
  query_part_2 <- ";\n  wdt:P225 ?taxaLabels. \n  ?structure wdt:P235 ?structure_id;\n  wdt:P233 ?structureSmiles;\n  p:P703 ?statement.\n  ?statement ps:P703 ?taxa;\n  prov:wasDerivedFrom ?ref.\n  ?ref pr:P248 ?references.\n  SERVICE <https://query-scholarly.wikidata.org/sparql> { \n ?references wdt:P1476 ?referencesLabels;\n  wdt:P356 ?references_ids;\n  wdt:P577 ?art_date.\n  }\n FILTER(((YEAR(?art_date)) >= "
  query_part_3 <- " ) && ((YEAR(?art_date)) <= "
  query_part_4 <- " ))\n  SERVICE wikibase:label { bd:serviceParam wikibase:language \"[AUTO_LANGUAGE], mul, en\". }\n}"

  lotus_path <- "data/source/libraries/lotus.csv.gz"
  message("Getting last LOTUS version")
  get_zenodo_file(
    doi = "10.5281/zenodo.5794106",
    pattern = "frozen_metadata.csv.gz",
    path = lotus_path
  )

  message("Loading LOTUS classified structures")
  lotus_raw <- fread_gz(
    file = lotus_path,
    select = c(
      "structure_inchikey",
      "structure_smiles_2D",
      "structure_exact_mass",
      "structure_xlogp",
      "structure_taxonomy_npclassifier_01pathway",
      "structure_taxonomy_npclassifier_02superclass",
      "structure_taxonomy_npclassifier_03class"
    )
  )
  structures_classified <- data.table::copy(lotus_raw)
  data.table::setnames(
    structures_classified,
    old = c(
      "structure_inchikey",
      "structure_taxonomy_npclassifier_01pathway",
      "structure_taxonomy_npclassifier_02superclass",
      "structure_taxonomy_npclassifier_03class"
    ),
    new = c(
      "structure_id",
      "chemical_pathway",
      "chemical_superclass",
      "chemical_class"
    )
  )
  structures_classified <- unique(structures_classified)

  message("Converting taxa to QIDs")
  qids <- stats::setNames(lapply(taxa, taxon_name_to_qid), taxa)

  message("Building queries")
  queries <- queries_progress(
    xs = qids,
    start = start,
    end = end,
    limit = limit,
    query_part_1 = query_part_1,
    query_part_2 = query_part_2,
    query_part_3 = query_part_3,
    query_part_4 = query_part_4
  )

  message("Querying Wikidata")
  results <- wiki_progress(xs = queries)

  message("Removing empty results")
  results <- Filter(function(res) nrow(res) > 0L, results)

  message("Cleaning tables and adding columns")
  tables <- tables_progress(
    xs = results,
    structures_classified = structures_classified
  )

  if (no_stereo) {
    tables <- lapply(tables, make_no_stereo)
  }

  if (filter_ms_conditions) {
    tables <- lapply(tables, make_chromatographiable)
  }

  message("Generating chemical hierarchies...")
  message("... for single taxa")
  hierarchies_simple <- hierarchies_progress(tables, comparison = comparison)
  message("... for grouped taxa")
  hierarchies_grouped <- hierarchies_grouped_progress(xs = tables)
  message("... combining")
  if (length(hierarchies_grouped) > 0L) {
    names(hierarchies_grouped) <- ifelse(
      !grepl("\\W+", names(hierarchies_grouped)),
      paste(names(hierarchies_grouped), "grouped", sep = "_"),
      names(hierarchies_grouped)
    )
  }
  hierarchies_grouped <- Filter(Negate(is.null), hierarchies_grouped)

  hierarchies <- c(hierarchies_simple, hierarchies_grouped)

  if (!is.null(comparison)) {
    message("Generating special comparison")
    special_list <- lapply(comparison, function(comp_nm) tables[[comp_nm]])
    spec_dt <- data.table::rbindlist(
      special_list,
      use.names = TRUE,
      fill = TRUE
    )
    has_multiword <- all(grepl("\\W+", comparison))
    spec_dt[, `:=`(
      best_candidate_1 = chemical_pathway,
      best_candidate_2 = chemical_superclass,
      best_candidate_3 = chemical_class,
      organism = if (has_multiword) taxaLabels else gsub(" .*", "", taxaLabels)
    )]
    spec_dt[, `:=`(sample = organism, species = organism)]
    drop_cols <- intersect(
      c("taxa", "taxaLabels", "references", "referencesLabels"),
      names(spec_dt)
    )
    if (length(drop_cols) > 0L) {
      spec_dt[, (drop_cols) := NULL]
    }
    hierarchies[["special"]] <- prepare_hierarchy(
      dataframe = unique(spec_dt),
      type = "literature"
    )
  }

  prepared_plots <- lapply(hierarchies, prepare_plot)
  plots <- lapply(prepared_plots, plot_histograms_litt, label = "")

  non_grouped_names <- names(hierarchies)[
    !grepl("_grouped", names(hierarchies))
  ]

  message("Generating treemaps")
  treemaps <- if (show_title) {
    treemaps_progress(
      xs = non_grouped_names,
      hierarchies = hierarchies,
      tables = tables
    )
  } else {
    treemaps_progress_no_title(
      xs = non_grouped_names,
      hierarchies = hierarchies
    )
  }

  message("Generating sunbursts")
  sunbursts <- if (show_title) {
    treemaps_progress(
      xs = non_grouped_names,
      type = "sunburst",
      hierarchies = hierarchies,
      tables = tables
    )
  } else {
    treemaps_progress_no_title(
      xs = non_grouped_names,
      type = "sunburst",
      hierarchies = hierarchies
    )
  }

  message("Filtering treemaps")
  treemaps <- treemaps[!grepl("ae$", names(treemaps))]

  message("Filtering sunbursts")
  sunbursts <- sunbursts[!grepl("ae$", names(sunbursts))]

  tables_out <- if (format_tables) {
    stats::setNames(
      lapply(names(tables), function(nm) format_gt(tables[[nm]], title = nm)),
      names(tables)
    )
  } else {
    tables
  }

  list(
    plots = plots,
    treemaps = treemaps,
    sunbursts = sunbursts,
    tables = tables_out
  )
}
