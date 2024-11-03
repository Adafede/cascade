start <- Sys.time()

pkgload::load_all()

source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/create_dir.R")
source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/get_file.R")
source(file = "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima/main/R/get_last_version_from_zenodo.R")

message("This program plots IDs")
message("Authors: \n", "AR")
message("Contributors: \n", "...")

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

message("Generating treemaps")
treemaps <-
  treemaps_progress_no_title(xs = names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))])

message("Generating sunbursts")
sunbursts <-
  treemaps_progress_no_title(xs = names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))], type = "sunburst")

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
