start <- Sys.time()

library(package = dplyr, quietly = TRUE)
library(package = gt, quietly = TRUE)
library(package = htmltools, quietly = TRUE)
library(package = microshades, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = RCurl, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = splitstackshape, quietly = TRUE)
library(package = tidyr, quietly = TRUE)
library(package = WikidataQueryServiceR, quietly = TRUE)

source(file = "R/colors.R")
source(file = "R/format_gt.R")
source(file = "R/plot_histograms.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_plot.R")
source(file = "R/log_debug.R")
source(file = "R/molinfo.R")

classified_path <-
  "~/Git/lotus-processor/data/processed/210715_frozen_metadata.csv.gz"

query_path_1 <- "src/sparql/get_review_table_part_1.sql"
query_path_2 <- "src/sparql/get_review_table_part_2.sql"
query_path_3 <- "src/sparql/get_review_table_part_3.sql"
query_path_4 <- "src/sparql/get_review_table_part_4.sql"

export_dir <- "data"

## As there is no better way than to manually assess if the QID
## really corresponds to what you want
qids <- c(
  # "Actinobacteria" = "Q26262282",
  # "Simaroubaceae" = "Q156679",
  "Swertia" = "Q163970",
  "Gentiana lutea" = "Q158572"
)

limit <- 1000
start_date <- 1900
end_date <- 2021

clean_xanthones <- TRUE

query_part_1 <- readr::read_file(query_path_1)
query_part_2 <- readr::read_file(query_path_2)
query_part_3 <- readr::read_file(query_path_3)
query_part_4 <- readr::read_file(query_path_4)

structures_classified <- readr::read_delim(
  file = classified_path,
  col_select = c(
    "structure_id" = "structure_inchikey",
    "chemical_pathway" = "structure_taxonomy_npclassifier_01pathway",
    "chemical_superclass" = "structure_taxonomy_npclassifier_02superclass",
    "chemical_class" = "structure_taxonomy_npclassifier_03class"
  )
) %>%
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
# ) %>%
#   dplyr::distinct()

queries <- character()
for (i in names(qids)) {
  queries[[i]] <- paste0(
    query_part_1,
    qids[i],
    query_part_2,
    start_date,
    query_part_3,
    end_date,
    query_part_4,
    paste("\nLIMIT", limit)
  )
}

results <- list()
for (i in names(queries)) {
  results[[i]] <-
    WikidataQueryServiceR::query_wikidata(sparql_query = queries[i])
}

tables <- list()
for (i in names(results)) {
  tables[[i]] <- results[[i]] %>%
    dplyr::left_join(structures_classified) %>%
    dplyr::mutate(structureImage = RCurl::curlEscape(structureSmiles)) %>%
    dplyr::relocate(structureImage, .after = structure) %>%
    dplyr::relocate(structureLabel, .before = structure) %>%
    dplyr::select(-references_ids, -structure_id, -structureSmiles) %>%
    splitstackshape::cSplit(
      c("taxa", "taxaLabels", "references", "referencesLabels"),
      sep = "|",
      direction = "long"
    ) %>%
    dplyr::group_by(structure) %>%
    tidyr::fill(c("taxa", "taxaLabels", "references", "referencesLabels"),
      .direction = "downup"
    ) %>%
    dplyr::group_by(chemical_class) %>%
    dplyr::add_count(sort = TRUE) %>%
    dplyr::select(-n) %>%
    dplyr::group_by(chemical_superclass) %>%
    dplyr::add_count(sort = TRUE) %>%
    dplyr::select(-n) %>%
    dplyr::group_by(chemical_pathway) %>%
    dplyr::add_count(sort = TRUE) %>%
    dplyr::select(-n) %>%
    dplyr::distinct()
}

subtables <- list()
for (i in names(tables)) {
  subtables[[i]] <- tables[[i]] %>%
    dplyr::filter(chemical_pathway == .[1, "chemical_pathway"]) %>%
    dplyr::group_by(chemical_class) %>%
    dplyr::add_count(sort = TRUE) %>%
    dplyr::select(-n) %>%
    dplyr::group_by(chemical_superclass) %>%
    dplyr::add_count(sort = TRUE) %>%
    dplyr::select(-n, -chemical_pathway) %>%
    dplyr::distinct()
}

prettyTables <- list()
for (i in names(tables)) {
  prettyTables[[i]] <-
    temp_gt_function(
      table = tables[[i]],
      title = i,
      subtitle = "All compounds"
    )
}

prettySubtables <- list()
for (i in names(subtables)) {
  prettySubtables[[i]] <-
    temp_gt_function(
      table = subtables[[i]],
      title = i,
      subtitle = tables[[i]][1, "chemical_pathway"]
    )
}

for (i in names(prettyTables)) {
  gt::gtsave(data = prettyTables[[i]], filename = file.path(export_dir, paste0(
    "prettyTable_",
    gsub(
      pattern = " ",
      replacement = "_",
      x = i
    ),
    ".html"
  )))
}

for (i in names(prettySubtables)) {
  gt::gtsave(data = prettySubtables[[i]], filename = file.path(export_dir, paste0(
    "prettySubtable_",
    gsub(
      pattern = " ",
      replacement = "_",
      x = i
    ),
    ".html"
  )))
}

dataframe_genus <- tables$`Swertia` |>
  mutate(
    best_candidate_1 = chemical_pathway,
    best_candidate_2 = chemical_superclass,
    best_candidate_3 = chemical_class,
    organism = taxaLabels,
    sample = taxaLabels,
    species = taxaLabels
  )

dataframe_species <- tables$`Gentiana lutea` |>
  mutate(
    best_candidate_1 = chemical_pathway,
    best_candidate_2 = chemical_superclass,
    best_candidate_3 = chemical_class,
    organism = taxaLabels,
    sample = taxaLabels,
    species = taxaLabels
  )

test_genus <- dataframe_genus |>
  prepare_hierarchy(type = "literature")

test_species <- dataframe_species |>
  prepare_hierarchy(type = "literature")

test_genus_2 <-
  prepare_plot(dataframe = test_genus)

test_species_2 <-
  prepare_plot(dataframe = test_species)

plot_histograms(
  dataframe = test_genus_2,
  label = "Based on literature"
)

plot_histograms(
  dataframe = test_species_2,
  label = "Based on literature"
)

plotly::plot_ly(
  data = test_genus |>
    filter(sample == "Swertia japonica"),
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "sunburst",
  branchvalues = "total"
) |>
  plotly::layout(colorway = sunburst_colors)

plotly::plot_ly(
  data = test_species,
  ids = ~ids,
  labels = ~labels,
  parents = ~parents,
  values = ~values,
  maxdepth = 3,
  type = "sunburst",
  branchvalues = "total"
) |>
  plotly::layout(colorway = sunburst_colors)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
