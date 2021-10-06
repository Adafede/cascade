library(package = dplyr, quietly = TRUE)
library(package = gt, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = RCurl, quietly = TRUE)
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

export_table_1_path <- "data/prettyTable.html"
export_subtable_1_path <- "data/prettySubtable.html"


query_part_1 <- read_file(query_path_1)
query_part_2 <- read_file(query_path_2)
query_part_3 <- read_file(query_path_3)
query_part_4 <- read_file(query_path_4)
qid <- "Q26262282"
start_date <- 1900
end_date <- 2021

myquery <- paste0(query_part_1,
                  qid,
                  query_part_2,
                  start_date,
                  query_part_3,
                  end_date,
                  query_part_4)

results <- query_wikidata(sparql_query = myquery)

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

table <- results %>%
  sample_n(100) %>%
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
       .direction = "downup") %>%
  group_by(chemical_pathway) %>%
  add_count(sort = TRUE) %>%
  select(-n) %>%
  distinct()

subtable_chemical <- table %>%
  filter(chemical_pathway == .[1, "chemical_pathway"]) %>%
  group_by(chemical_superclass) %>%
  add_count(sort = TRUE) %>%
  select(-n, -chemical_pathway) %>%
  distinct()

prettyTable <- temp_gt_function(table = table)

prettySubtable <- temp_gt_function(table = subtable_chemical)

gtsave(data = prettyTable, filename = export_table_1_path)
gtsave(data = prettySubtable, filename = export_subtable_1_path)
