#' Generate IDs
#'
#' @export
#'
#' @include hierarchies_progress.R
#' @include hierarchies_grouped_progress.R
#' @include make_chromatographiable.R
#' @include make_no_stereo.R
#' @include plot_histograms.R
#' @include prepare_hierarchy.R
#' @include prepare_plot.R
#' @include queries_progress.R
#' @include tables_progress.R
#' @include taxon_name_to_qid.R
#' @include treemaps_progress.R
#' @include wiki_progress.R
#'
#' @param taxa Taxa
#' @param comparison Comparison
#' @param no_stereo No stereo
#' @param filter_ms_conditions Filter MS conditions
#' @param start Start
#' @param end End
#' @param limit Limit
#' @param create_dir_f Create dir function
#' @param get_file_f Get file function
#' @param get_last_version_from_zenodo_f Get last version from Zenodo function
#'
#' @return IDs
#'
#' @examples
#' \dontrun{
#' generate_ids()
#' }
generate_ids <- function(taxa = c("Swertia", "Kopsia", "Ginkgo"),
                         comparison = c("Swertia", "Kopsia"),
                         no_stereo = TRUE,
                         filter_ms_conditions = TRUE,
                         start = "0",
                         end = "9999",
                         limit = "1000000",
                         create_dir_f = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/create_dir.R",
                         get_file_f = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/get_file.R",
                         get_last_version_from_zenodo_f = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/get_last_version_from_zenodo.R") {
  source(file = create_dir_f)
  source(file = get_file_f)
  source(file = get_last_version_from_zenodo_f)
  query_part_1 <- "SELECT ?structure ?structureLabel ?structure_id ?structureSmiles (GROUP_CONCAT(?taxon_name; SEPARATOR = \"|\") AS ?taxaLabels) (GROUP_CONCAT(?taxon; SEPARATOR = \"|\") AS ?taxa) (GROUP_CONCAT(?art_title; SEPARATOR = \"|\") AS ?referencesLabels) (GROUP_CONCAT(?art_doi; SEPARATOR = \"|\") AS ?references_ids) (GROUP_CONCAT(?art; SEPARATOR = \"|\") AS ?references) WHERE {\n  ?taxon (wdt:P171*) wd:"
  query_part_2 <- ";\n  wdt:P225 ?taxon_name.\n  ?structure wdt:P235 ?structure_id;\n  wdt:P233 ?structureSmiles;\n  p:P703 ?statement.\n  ?statement ps:P703 ?taxon;\n  prov:wasDerivedFrom ?ref.\n  ?ref pr:P248 ?art.\n  ?art wdt:P1476 ?art_title;\n  wdt:P356 ?art_doi;\n  wdt:P577 ?art_date.\n  FILTER(((YEAR(?art_date)) >= "
  query_part_3 <- " ) && ((YEAR(?art_date)) <= "
  query_part_4 <- " ))\n  SERVICE wikibase:label { bd:serviceParam wikibase:language \"[AUTO_LANGUAGE],en\". }\n}\nGROUP BY ?structure ?structure_id ?structureLabel ?structureSmiles"

  message("Getting last LOTUS version")
  get_last_version_from_zenodo(
    doi = "10.5281/zenodo.5794106",
    pattern = "frozen_metadata.csv.gz",
    "data/source/libraries/lotus.csv.gz"
  )

  message("Loading LOTUS classified structures")
  structures_classified <- tidytable::fread(
    file = "data/source/libraries/lotus.csv.gz",
    select = c(
      "structure_inchikey",
      "structure_smiles_2D",
      "structure_exact_mass",
      "structure_xlogp",
      "structure_taxonomy_npclassifier_01pathway",
      "structure_taxonomy_npclassifier_02superclass",
      "structure_taxonomy_npclassifier_03class"
    )
  ) |>
    tidytable::rename(
      "structure_id" = "structure_inchikey",
      "chemical_pathway" = "structure_taxonomy_npclassifier_01pathway",
      "chemical_superclass" = "structure_taxonomy_npclassifier_02superclass",
      "chemical_class" = "structure_taxonomy_npclassifier_03class"
    ) |>
    tidytable::distinct()

  message("Converting taxa to QIDs")
  qids <- taxa |>
    lapply(taxon_name_to_qid)
  names(qids) <- taxa
  comparison <- c("Swertia", "Kopsia")

  genera <-
    names(qids)[!grepl(
      pattern = " ",
      x = names(qids),
      fixed = TRUE
    )]

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
  results <- purrr::keep(results, ~ nrow(.) > 0)

  message("Cleaning tables and adding columns")
  tables <- tables_progress(xs = results, structures_classified = structures_classified)

  if (no_stereo) {
    tables <- lapply(tables, make_no_stereo)
  }

  if (filter_ms_conditions) {
    tables <- lapply(tables, make_chromatographiable)
  }

  # tables$SwertiaExp <-
  #   readr::read_tsv("data/paper/prettyTable.tsv") |>
  #   dplyr::select(
  #     structure = Structure,
  #     chemical_pathway = `Chemical Pathway`,
  #     chemical_superclass = `Chemical Superclass`,
  #     chemical_class = `Chemical Class`,
  #     n = `Peak Area [%]`
  #   ) |>
  #   dplyr::mutate(
  #     structureLabel = "bla",
  #     structureImage = "bla",
  #     taxaLabels = "SwertiaExp",
  #     taxa = "SwertiaExp",
  #     referencesLabels = "bla",
  #     references = "bla",
  #     structure_exact_mass = 0,
  #     structure_xlop = 0,
  #   ) |>
  #   dplyr::rowwise() |>
  #   dplyr::slice(rep(1:dplyr::n(), each = n * 100)) |>
  #   dplyr::filter(!is.na(structure))

  message("Generating chemical hierarchies...")
  message("... for single taxa")
  hierarchies_simple <-
    hierarchies_progress(tables)
  message("... for grouped taxa")
  hierarchies_grouped <- hierarchies_grouped_progress(xs = tables)
  message("... combining")
  names(hierarchies_grouped) <- ifelse(
    test = !grepl(pattern = "\\W+", x = names(hierarchies_grouped)),
    yes = paste(names(hierarchies_grouped), "grouped", sep = "_"),
    no = names(hierarchies_grouped)
  )
  hierarchies_grouped <- purrr::compact(hierarchies_grouped)

  hierarchies <- append(hierarchies_simple, hierarchies_grouped)

  if (!is.null(comparison)) {
    message("Generating special comparison")
    special <- lapply(
      X = seq_along(comparison),
      FUN = function(x) {
        tables[[comparison[[x]]]]
      }
    )
    hierarchies[["special"]] <- prepare_hierarchy(
      dataframe = dplyr::bind_rows(special) |>
        dplyr::rowwise() |>
        dplyr::mutate(
          best_candidate_1 = chemical_pathway,
          best_candidate_2 = chemical_superclass,
          best_candidate_3 = chemical_class,
          organism = ifelse(
            test = all(grepl(pattern = "\\W+", x = comparison)),
            yes = taxaLabels,
            no = gsub(
              pattern = " .*",
              replacement = "",
              x = taxaLabels
            )
          )
        ) |>
        dplyr::mutate(sample = organism, species = organism) |>
        dplyr::select(-taxa, -taxaLabels, -references, -referencesLabels) |>
        dplyr::distinct(),
      type = "literature"
    )
  }

  prepared_plots <- hierarchies |>
    lapply(prepare_plot)

  plots <- prepared_plots |>
    lapply(plot_histograms_litt, label = "")

  message("Generating treemaps")
  treemaps <-
    treemaps_progress_no_title(xs = names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))], hierarchies = hierarchies)

  message("Generating sunbursts")
  sunbursts <-
    treemaps_progress_no_title(
      xs = names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))],
      type = "sunburst",
      hierarchies = hierarchies
    )

  message("Filtering treemaps")
  treemaps <-
    within(treemaps, rm(list = names(treemaps)[grepl(pattern = "ae$", x = names(treemaps))]))

  message("Filtering sunbursts")
  sunbursts <-
    within(sunbursts, rm(list = names(sunbursts)[grepl(pattern = "ae$", x = names(sunbursts))]))

  return(list(
    plots = plots,
    treemaps = treemaps,
    sunbursts = sunbursts
  ))
}
