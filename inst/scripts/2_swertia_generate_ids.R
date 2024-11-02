start <- Sys.time()

source(file = "R/check_and_load_packages.R")
source(file = "R/check_export_dir.R")
source(file = "R/colors.R")
source(file = "R/format_gt.R")
source(file = "R/hierarchies_progress.R")
source(file = "R/hierarchies_grouped_progress.R")
source(file = "R/histograms_progress.R")
source(file = "R/make_2D.R")
source(file = "R/make_chromatographiable.R")
source(file = "R/molinfo.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/plot_histograms.R")
source(file = "R/prehistograms_progress.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_plot.R")
source(file = "R/prettyTables_progress.R")
source(file = "R/queries_progress.R")
source(file = "R/save_histograms_progress.R")
source(file = "R/save_prettySubtables_progress.R")
source(file = "R/save_prettyTables_progress.R")
source(file = "R/save_treemaps_progress.R")
source(file = "R/subtables_progress.R")
source(file = "R/tables_progress.R")
source(file = "R/treemaps_progress.R")
source(file = "R/wiki_progress.R")
source(file = "R/cascade-package.R")
source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/create_dir.R")
source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/get_file.R")
source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/get_last_version_from_zenodo.R")

paths <- parse_yaml_paths()
get_last_version_from_zenodo(
  doi = paths$url$lotus$doi,
  pattern = paths$urls$lotus$pattern,
  path = paths$data$source$libraries$lotus
)

exports <-
  list(
    paths$data$path,
    paths$data$histograms$path,
    paths$data$sunbursts$path,
    paths$data$tables$path,
    paths$data$treemaps$path
  )

taxa <- c("Gentiana", "Kopsia", "Ginkgo")
dimensionality <- 2
c18 <- TRUE

taxon_name_to_qid <- function(taxon_name) {
  WikidataQueryServiceR::query_wikidata(
    sparql_query = paste0(
      "SELECT ?search ?item WHERE {
    SERVICE wikibase:mwapi {
      bd:serviceParam wikibase:endpoint \"www.wikidata.org\";
                      wikibase:api \"EntitySearch\";
                      mwapi:search \"",
      taxon_name,
      "                 \";
                      mwapi:language \"mul\".
    ?item wikibase:apiOutputItem mwapi:item.
    ?num wikibase:apiOrdinal true.
    }
    ?item (wdt:P225) ?search.
    FILTER (?num = 0)
    }"
    )
  ) |>
    tidytable::pull(item) |>
    gsub(pattern = "http://www.wikidata.org/entity/", replacement = "")
}

#' As there is no better way than to manually assess if the QID
#' really corresponds to what you want
qids <- taxa |>
  lapply(taxon_name_to_qid)
names(qids) <- taxa
comparison <- taxa <- c("Gentiana", "Alstonia")

genera <-
  names(qids)[!grepl(
    pattern = " ",
    x = names(qids),
    fixed = TRUE
  )]

query_part_1 <- readr::read_file(paths$inst$scripts$sparql$review_1)
query_part_2 <- readr::read_file(paths$inst$scripts$sparql$review_2)
query_part_3 <- readr::read_file(paths$inst$scripts$sparql$review_3)
query_part_4 <- readr::read_file(paths$inst$scripts$sparql$review_4)

message("Loading LOTUS classified structures")
structures_classified <- readr::read_delim(
  file = paths$data$source$libraries$lotus,
  col_select = c(
    "structure_id" = "structure_inchikey",
    "structure_smiles_2D",
    "structure_exact_mass",
    "structure_xlogp",
    "chemical_pathway" = "structure_taxonomy_npclassifier_01pathway",
    "chemical_superclass" = "structure_taxonomy_npclassifier_02superclass",
    "chemical_class" = "structure_taxonomy_npclassifier_03class"
  )
) |>
  dplyr::distinct()

# organisms_classified <- readr::read_delim(
#   file = classified_path,
#   col_select = c(
#     "organism" = "organism_wikidata",
#     "taxonLabel" = "organism_name",
#     "taxonId" = "organism_taxonomy_ottid",
#     "taxon_01domain" = "organism_taxonomy_01domain",
#     "taxon_02kingdom" = "organism_taxonomy_02kingdom",
#     "taxon_03phylum" = "organism_taxonomy_03phylum",
#     "taxon_04class" = "organism_taxonomy_04class",
#     "taxon_05order" = "organism_taxonomy_05order",
#     "taxon_06family" = "organism_taxonomy_06family",
#     "taxon_07tribe" = "organism_taxonomy_07tribe",
#     "taxon_08genus" = "organism_taxonomy_08genus",
#     "taxon_09species" = "organism_taxonomy_09species",
#     "taxon_10varietas" = "organism_taxonomy_10varietas"
#   )
# ) |>
#   dplyr::distinct()

message("Building queries")
queries <- queries_progress(xs = qids)

message("Querying Wikidata")
results <- wiki_progress(xs = queries)

message("Removing empty results")
results <- purrr::keep(results, ~ nrow(.) > 0)

message("Cleaning tables and adding columns")
tables <- tables_progress(xs = results)

if (dimensionality == 2) {
  tables <- lapply(tables, make_2D)
}

if (c18 == TRUE) {
  tables <- lapply(tables, make_chromatographiable)
}

tables$SwertiaExp <-
  readr::read_tsv("data/paper/prettyTable.tsv") |>
  dplyr::select(
    structure = Structure,
    chemical_pathway = `Chemical Pathway`,
    chemical_superclass = `Chemical Superclass`,
    chemical_class = `Chemical Class`,
    n = `Peak Area [%]`
  ) |>
  dplyr::mutate(
    structureLabel = "bla",
    structureImage = "bla",
    taxaLabels = "SwertiaExp",
    taxa = "SwertiaExp",
    referencesLabels = "bla",
    references = "bla",
    structure_exact_mass = 0,
    structure_xlop = 0,
  ) |>
  dplyr::rowwise() |>
  dplyr::slice(rep(1:dplyr::n(), each = n * 100)) |>
  dplyr::filter(!is.na(structure))

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

#' Title
#'
#' @param xs
#' @param type
#'
#' @return
#' @export
#'
#' @examples
treemaps_progress <- function(xs, type = "treemap") {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = setNames(object = xs, nm = xs),
    FUN = function(x) {
      p()
      if (x != "special") {
        plotly::plot_ly(
          data = hierarchies[[x]],
          ids = ~ids,
          labels = ~labels,
          parents = ~parents,
          values = ~values,
          maxdepth = 3,
          type = type,
          branchvalues = "total",
          textinfo = "label+percent value+percent parent+percent root"
        ) |>
          plotly::layout(
            colorway = microshades_colors,
            title = paste(x, "(", nrow(
              tables[[x]] |> dplyr::distinct(structure)
            ), ")"),
            margin = list(t = 40)
          )
      } else {
        plotly::plot_ly() |>
          plotly::add_trace(
            data = hierarchies[[unique(hierarchies[[x]]$species)[1]]],
            ids = ~ids,
            labels = ~labels,
            parents = ~parents,
            values = ~values,
            maxdepth = 3,
            type = type,
            branchvalues = "total",
            textinfo = "label+percent value+percent parent+percent root",
            domain = list(row = 0, column = 0)
          ) |>
          plotly::add_trace(
            data = hierarchies[[unique(hierarchies[[x]]$species)[2]]],
            ids = ~ids,
            labels = ~labels,
            parents = ~parents,
            values = ~values,
            maxdepth = 3,
            type = type,
            branchvalues = "total",
            textinfo = "label+percent value+percent parent+percent root",
            domain = list(row = 0, column = 1)
          ) |>
          plotly::layout(
            # title = paste(
            #   "Comparative analysis",
            #   "\n",
            #   unique(hierarchies[[x]]$species)[1],
            #   "(",
            #   nrow(tables[[unique(hierarchies[[x]]$species)[1]]] |> dplyr::distinct(structure)),
            #   ")",
            #   "                                 ",
            #   unique(hierarchies[[x]]$species)[2],
            #   "(",
            #   nrow(tables[[unique(hierarchies[[x]]$species)[2]]] |> dplyr::distinct(structure)),
            #   ")"
            # ),
            grid = list(rows = 1, columns = 2),
            colorway = microshades_colors,
            margin = list(t = 40)
          )
      }
    }
  )
}
message("Generating treemaps")
treemaps <-
  treemaps_progress(xs = names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))])

message("Generating sunbursts")
sunbursts <-
  treemaps_progress(xs = names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))], type = "sunburst")

message("Filtering treemaps")
treemaps <-
  within(treemaps, rm(list = names(treemaps)[grepl(pattern = "ae$", x = names(treemaps))]))

message("Filtering sunbursts")
sunbursts <-
  within(sunbursts, rm(list = names(sunbursts)[grepl(pattern = "ae$", x = names(sunbursts))]))

# plotly::save_image(
#   p = treemaps$special,
#   file = "data/paper/cascade-8.pdf",
#   width = 900,
#   height = 900
# )

end <- Sys.time()

message("Script finished in ", format(end - start))
