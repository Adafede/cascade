library(package = dplyr, quietly = TRUE)
library(package = gt, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = RCurl, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = splitstackshape, quietly = TRUE)
library(package = tidyr, quietly = TRUE)
library(package = WikidataQueryServiceR, quietly = TRUE)

source(file = "R/molinfo.R")
source(file = "R/format_gt.R")

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
  "Actinobacteria" = "Q26262282",
  "Simaroubaceae" = "Q156679",
  "Swertia" = "Q163970",
  "Gentiana lutea" = "Q158572"
)

start_date <- 1900
end_date <- 2021

query_part_1 <- read_file(query_path_1)
query_part_2 <- read_file(query_path_2)
query_part_3 <- read_file(query_path_3)
query_part_4 <- read_file(query_path_4)

classified <- read_delim(
  file = classified_path,
  col_select = c(
    "structure_id" = "structure_inchikey",
    "chemical_pathway" = "structure_taxonomy_npclassifier_01pathway",
    "chemical_superclass" = "structure_taxonomy_npclassifier_02superclass",
    "chemical_class" = "structure_taxonomy_npclassifier_03class"
  )
) %>%
  distinct()

queries <- character()
for (i in names(qids)) {
  queries[[i]] <- paste0(
    query_part_1,
    qids[i],
    query_part_2,
    start_date,
    query_part_3,
    end_date,
    query_part_4
  )
}

results <- list()
for (i in names(queries)) {
  results[[i]] <- query_wikidata(sparql_query = queries[i])
}

tables <- list()
for (i in names(results)) {
  tables[[i]] <- results[[i]] %>%
    left_join(classified) %>%
    mutate(structureImage = curlEscape(structureSmiles)) %>%
    relocate(structureImage, .after = structure) %>%
    relocate(structureLabel, .before = structure) %>%
    select(-references_ids, -structure_id, -structureSmiles) %>%
    cSplit(
      c("taxa", "taxaLabels", "references", "referencesLabels"),
      sep = "|",
      direction = "long"
    ) %>%
    group_by(structure) %>%
    fill(c("taxa", "taxaLabels", "references", "referencesLabels"),
      .direction = "downup"
    ) %>%
    group_by(chemical_class) %>%
    add_count(sort = TRUE) %>%
    select(-n) %>%
    group_by(chemical_superclass) %>%
    add_count(sort = TRUE) %>%
    select(-n) %>%
    group_by(chemical_pathway) %>%
    add_count(sort = TRUE) %>%
    select(-n) %>%
    distinct()
}

subtables <- list()
for (i in names(tables)) {
  subtables[[i]] <- tables[[i]] %>%
    filter(chemical_pathway == .[1, "chemical_pathway"]) %>%
    group_by(chemical_class) %>%
    add_count(sort = TRUE) %>%
    select(-n) %>%
    group_by(chemical_superclass) %>%
    add_count(sort = TRUE) %>%
    select(-n, -chemical_pathway) %>%
    distinct()
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
  gtsave(data = prettyTables[[i]], filename = file.path(export_dir, paste0(
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
  gtsave(data = prettySubtables[[i]], filename = file.path(export_dir, paste0(
    "prettySubtable_",
    gsub(
      pattern = " ",
      replacement = "_",
      x = i
    ),
    ".html"
  )))
}
