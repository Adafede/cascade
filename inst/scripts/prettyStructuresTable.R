start <- Sys.time()

library(package = dplyr, quietly = TRUE)
library(package = ggplot2, quietly = TRUE)
library(package = gt, quietly = TRUE)
library(package = htmltools, quietly = TRUE)
library(package = microshades, quietly = TRUE)
library(package = plotly, quietly = TRUE)
library(package = purrr, quietly = TRUE)
library(package = RCurl, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = splitstackshape, quietly = TRUE)
library(package = tidyr, quietly = TRUE)
library(package = WikidataQueryServiceR, quietly = TRUE)

source(file = "R/check_export_dir.R")
source(file = "R/colors.R")
source(file = "R/format_gt.R")
source(file = "R/plot_histograms.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_plot.R")
source(file = "R/log_debug.R")
source(file = "R/molinfo.R")

classified_path <-
  "~/Git/lotus-processor/data/processed/211220_frozen_metadata.csv.gz"

query_path_1 <- "inst/scripts/sparql/get_review_table_part_1.sql"
query_path_2 <- "inst/scripts/sparql/get_review_table_part_2.sql"
query_path_3 <- "inst/scripts/sparql/get_review_table_part_3.sql"
query_path_4 <- "inst/scripts/sparql/get_review_table_part_4.sql"

export_dir <- "data"
export_dir_histograms <- file.path(export_dir, "histograms")
export_dir_sunbursts <- file.path(export_dir, "sunbursts")
export_dir_tables <- file.path(export_dir, "tables")
export_dir_treemaps <- file.path(export_dir, "treemaps")

exports <-
  list(
    export_dir,
    export_dir_histograms,
    export_dir_sunbursts,
    export_dir_tables,
    export_dir_treemaps
  )

#' As there is no better way than to manually assess if the QID
#' really corresponds to what you want
qids <- c(
  # "Actinobacteria" = "Q26262282"
  # "Asteraceae" = "Q25400"
  # "Simaroubaceae" = "Q156679"
  # "Gentianaceae" = "Q157216"
  # "Picrasma" = "Q135638",
  # "Picrasma quassioides" = "Q855778",
  "Swertia" = "Q163970",
  "Swertia chirayita" = "Q21318003",
  "Gentiana" = "Q144682",
  "Gentiana lutea" = "Q158572",
  "Quassia" = "Q1947702",
  "Quassia amara" = "Q135389"
  # "Dendrobium" = "Q133778",
  # "Dendrobium chrysanthum" = "Q5223343",
  # "Dendrobium fimbriatum" = "Q7990065",
  # "Trichoderma" = "Q135322",
  # "Trichoderma yunnanense" = "Q108442404",
  # "Papiliotrema" = "Q7132982",
  # "Papiliotrema rajasthanensis" = "Q27866418"
)

comparison <- c("Gentiana", "Swertia")
# comparison <- c("Dendrobium", "Trichoderma")

limit <- 50000
start_date <- 1900
end_date <- 2022

genera <-
  names(qids)[!grepl(
    pattern = " ",
    x = names(qids),
    fixed = TRUE
  )]

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
# ) |>
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
  if (nrow(results[[i]] != 0)) {
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
  } else {
    tables[[i]] <- data.frame() |>
      dplyr::mutate(
        structureLabel = NA,
        structure = NA,
        structureImage = NA,
        taxaLabels = NA,
        taxa = NA,
        referenceLabels = NA,
        references = NA,
        chemical_pathway = NA,
        chemical_superclass = NA,
        chemical_class = NA
      )
  }
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

hierarchies <- list()
for (i in names(tables)) {
  if (nrow(tables[[i]]) != 0) {
    hierarchies[[i]] <- prepare_hierarchy(
      dataframe = tables[[i]] |>
        dplyr::mutate(
          best_candidate_1 = chemical_pathway,
          best_candidate_2 = chemical_superclass,
          best_candidate_3 = chemical_class,
          organism = ifelse(
            test = grepl(pattern = "\\W+", x = i),
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
  } else {
    hierarchies[[i]] <- data.frame() |>
      dplyr::mutate(
        parents = NA,
        ids = NA,
        labels = NA,
        sample = NA,
        species = NA,
        values = 0
      ) |>
      dplyr::mutate_if(is.logical, as.character)
  }
}
for (i in names(tables)) {
  if (nrow(tables[[i]]) != 0) {
    if (!grepl(pattern = "\\W+", x = i)) {
      hierarchies[[paste0(i, "_grouped")]] <- prepare_hierarchy(
        dataframe = tables[[i]] |>
          dplyr::mutate(
            best_candidate_1 = chemical_pathway,
            best_candidate_2 = chemical_superclass,
            best_candidate_3 = chemical_class,
            organism = taxaLabels
          ) |>
          dplyr::mutate(sample = organism, species = organism),
        type = "literature"
      )
    }
  } else {
    hierarchies[[i]] <- data.frame() |>
      dplyr::mutate(
        parents = NA,
        ids = NA,
        labels = NA,
        sample = NA,
        species = NA,
        values = 0
      ) |>
      dplyr::mutate_if(is.logical, as.character)
  }
}

prehistograms <- list()
for (i in names(hierarchies)) {
  if (nrow(hierarchies[[i]]) != 0) {
    prehistograms[[i]] <-
      prepare_plot(dataframe = hierarchies[[i]])
  } else {
    prehistograms[[i]] <- data.frame() |>
      dplyr::mutate(
        parents = NA,
        ids = NA,
        labels = NA,
        sample = NA,
        species = NA,
        values = NA,
        group = NA,
        subgroup = NA,
        tot = NA,
        color = NA,
        relative = NA
      ) |>
      dplyr::rowwise()
  }
}

histograms <- list()
for (i in names(prehistograms)[!grepl(pattern = "\\W+", x = names(prehistograms))]) {
  histograms[[i]] <- plot_histograms(
    dataframe = prehistograms[[i]],
    label = "Organism"
  )
}

#' Special requests
special <- list()
for (i in (seq_along(1:length(comparison)))) {
  special[[i]] <- tables[[comparison[[i]]]]
}
hierarchies[["special"]] <- prepare_hierarchy(
  dataframe = dplyr::bind_rows(special) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      best_candidate_1 = chemical_pathway,
      best_candidate_2 = chemical_superclass,
      best_candidate_3 = chemical_class,
      organism = ifelse(
        test = grepl(pattern = "\\W+", x = i),
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
prehistograms[["special"]] <-
  prepare_plot(dataframe = hierarchies[["special"]])
histograms[["special"]] <-
  plot_histograms(
    dataframe = prehistograms[["special"]],
    label = "Organism"
  )

treemaps <- list()
for (i in names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))]) {
  if (i != "special") {
    treemaps[[i]] <- plotly::plot_ly(
      data = hierarchies[[i]] |>
        filter(sample == i),
      ids = ~ids,
      labels = ~labels,
      parents = ~parents,
      values = ~values,
      maxdepth = 3,
      type = "treemap",
      branchvalues = "total",
      textinfo = "label+value+percent parent+percent root"
    ) |>
      plotly::layout(
        colorway = sunburst_colors,
        title = paste(i, "(", nrow(
          tables[[i]] |> dplyr::distinct(structure)
        ), ")"),
        margin = list(t = 40)
      )
  } else {
    treemaps[[i]] <- plotly::plot_ly() |>
      add_trace(
        data = hierarchies[[unique(hierarchies[[i]]$species)[1]]] |>
          filter(sample == unique(hierarchies[[i]]$species)[1] &
            !is.na(labels)),
        ids = ~ids,
        labels = ~labels,
        parents = ~parents,
        values = ~values,
        maxdepth = 3,
        type = "treemap",
        branchvalues = "total",
        textinfo = "label+value+percent parent+percent root",
        domain = list(row = 0, column = 0)
      ) |>
      add_trace(
        data = hierarchies[[unique(hierarchies[[i]]$species)[2]]] |>
          filter(sample == unique(hierarchies[[i]]$species)[2] &
            !is.na(labels)),
        ids = ~ids,
        labels = ~labels,
        parents = ~parents,
        values = ~values,
        maxdepth = 3,
        type = "treemap",
        branchvalues = "total",
        textinfo = "label+value+percent parent+percent root",
        domain = list(row = 0, column = 1)
      ) |>
      plotly::layout(
        title = paste(
          "Comparative analysis",
          "\n",
          unique(hierarchies[[i]]$species)[1],
          "(",
          nrow(tables[[unique(hierarchies[[i]]$species)[1]]] |> dplyr::distinct(structure)),
          ")",
          "                                 ",
          unique(hierarchies[[i]]$species)[2],
          "(",
          nrow(tables[[unique(hierarchies[[i]]$species)[2]]] |> dplyr::distinct(structure)),
          ")"
        ),
        grid = list(rows = 1, columns = 2),
        colorway = sunburst_colors,
        margin = list(t = 40)
      )
  }
}

sunbursts <- list()
for (i in names(hierarchies)[!grepl(pattern = "_grouped", x = names(hierarchies))]) {
  if (i != "special") {
    sunbursts[[i]] <- plotly::plot_ly(
      data = hierarchies[[i]] |>
        filter(sample == i),
      ids = ~ids,
      labels = ~labels,
      parents = ~parents,
      values = ~values,
      maxdepth = 3,
      type = "sunburst",
      branchvalues = "total"
    ) |>
      plotly::layout(
        colorway = sunburst_colors,
        title = paste(i, "(", nrow(
          tables[[i]] |> dplyr::distinct(structure)
        ), ")"),
        margin = list(t = 40)
      )
  } else {
    sunbursts[[i]] <- plotly::plot_ly() |>
      add_trace(
        data = hierarchies[[unique(hierarchies[[i]]$species)[1]]] |>
          filter(sample == unique(hierarchies[[i]]$species)[1] &
            !is.na(labels)),
        ids = ~ids,
        labels = ~labels,
        parents = ~parents,
        values = ~values,
        maxdepth = 3,
        type = "sunburst",
        branchvalues = "total",
        domain = list(row = 0, column = 0)
      ) |>
      add_trace(
        data = hierarchies[[unique(hierarchies[[i]]$species)[2]]] |>
          filter(sample == unique(hierarchies[[i]]$species)[2] &
            !is.na(labels)),
        ids = ~ids,
        labels = ~labels,
        parents = ~parents,
        values = ~values,
        maxdepth = 3,
        type = "sunburst",
        branchvalues = "total",
        domain = list(row = 0, column = 1)
      ) |>
      plotly::layout(
        title = paste(
          "Comparative analysis",
          "\n",
          unique(hierarchies[[i]]$species)[1],
          "(",
          nrow(tables[[unique(hierarchies[[i]]$species)[1]]] |> dplyr::distinct(structure)),
          ")",
          "                                 ",
          unique(hierarchies[[i]]$species)[2],
          "(",
          nrow(tables[[unique(hierarchies[[i]]$species)[2]]] |> dplyr::distinct(structure)),
          ")"
        ),
        grid = list(rows = 1, columns = 2),
        colorway = sunburst_colors,
        margin = list(t = 40)
      )
  }
}

lapply(X = exports, FUN = check_export_dir)

for (i in names(prettyTables)) {
  gt::gtsave(
    data = prettyTables[[i]],
    filename = file.path(export_dir_tables, paste0(
      "prettyTable_",
      gsub(
        pattern = " ",
        replacement = "_",
        x = i
      ),
      ".html"
    ))
  )
}

for (i in names(prettySubtables)) {
  gt::gtsave(
    data = prettySubtables[[i]],
    filename = file.path(
      export_dir_tables,
      paste0(
        "prettySubtable_",
        gsub(
          pattern = " ",
          replacement = "_",
          x = i
        ),
        ".html"
      )
    )
  )
}

dimensions <- list()
for (i in names(prehistograms)) {
  dimensions[[i]] <-
    nrow(prehistograms[[i]] |> dplyr::ungroup() |> dplyr::distinct(ids))
}

size <- list()
for (i in names(prehistograms)) {
  size[[i]] <-
    nrow(prehistograms[[i]] |> dplyr::ungroup() |> dplyr::distinct(species))
}

for (i in names(histograms)) {
  ggplot2::ggsave(
    plot = histograms[[i]],
    filename = file.path(
      export_dir_histograms,
      paste0(
        "histogram_",
        gsub(
          pattern = " ",
          replacement = "_",
          x = i
        ),
        ".pdf"
      )
    ),
    width = 16 * max((dimensions[[i]] / 30), 0.5),
    height = 9 * min(max((size[[i]] / 100), 1), 10) * max((dimensions[[i]] / 30), 0.5)
  )
}

for (i in names(sunbursts)) {
  plotly::save_image(
    p = sunbursts[[i]],
    file = file.path(
      export_dir_sunbursts,
      paste0(
        "sunburst_",
        gsub(
          pattern = " ",
          replacement = "_",
          x = i
        ),
        ".pdf"
      )
    ),
    width = 900,
    height = 900
  )
}

for (i in names(treemaps)) {
  plotly::save_image(
    p = treemaps[[i]],
    file = file.path(
      export_dir_treemaps,
      paste0(
        "treemap_",
        gsub(
          pattern = " ",
          replacement = "_",
          x = i
        ),
        ".pdf"
      )
    ),
    width = 900,
    height = 900
  )
}

end <- Sys.time()

log_debug("Script finished in", format(end - start))
