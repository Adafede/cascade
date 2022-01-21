source(file = "r/log_debug.R")
log_debug("This script plots an alternative magic tree.")

start <- Sys.time()

library(package = data.table, quietly = TRUE)
library(package = dplyr, quietly = TRUE)
library(package = forcats, quietly = TRUE)
library(package = ggplot2, quietly = TRUE)
library(package = ggtree, quietly = TRUE)
library(package = ggtreeExtra, quietly = TRUE)
library(package = ggstar, quietly = TRUE)
library(package = ggnewscale, quietly = TRUE)
library(package = microshades, quietly = TRUE)
library(package = readr, quietly = TRUE)
library(package = rotl, quietly = TRUE)
library(package = splitstackshape, quietly = TRUE)
library(package = tidyr, quietly = TRUE)
source(file = "r/colors.R")
source(file = "r/prepare_hierarchy.R")
source(file = "r/prepare_plot.R")

classified_path <-
  "~/Git/lotus-processor/data/processed/211220_frozen_metadata.csv.gz"

n_min <- 50

pairs_metadata <- readr::read_delim(file = classified_path) %>%
  data.table::data.table()

families_restricted <- pairs_metadata |>
  dplyr::filter(!is.na(organism_taxonomy_06family)) |>
  dplyr::group_by(organism_taxonomy_06family) |>
  dplyr::add_count() |>
  dplyr::ungroup() |>
  dplyr::filter(n >= n_min) |>
  dplyr::distinct(organism_taxonomy_06family)

families_matched_restricted <-
  rotl::tnrs_match_names(
    names = families_restricted$organism_taxonomy_06family,
    do_approximate_matching = FALSE
  )

ott_in_tree <-
  rotl::ott_id(families_matched_restricted)[rotl::is_in_tree(rotl::ott_id(families_matched_restricted))]

families_restricted <- families_restricted |>
  dplyr::filter(organism_taxonomy_06family %in% names(ott_in_tree))

families_matched_restricted <-
  rotl::tnrs_match_names(
    names = families_restricted$organism_taxonomy_06family,
    do_approximate_matching = FALSE
  )

families_matched_restricted <- families_matched_restricted |>
  dplyr::mutate(key = paste(
    gsub(
      x = unique_name,
      pattern = " ",
      replacement = "_",
      fixed = TRUE
    ),
    paste0("ott", ott_id),
    sep = "_"
  ))

tr_restricted <- rotl::tol_induced_subtree(ott_ids = ott_in_tree)

specific_classes <- pairs_metadata |>
  splitstackshape::cSplit(
    splitCols = colnames(pairs_metadata)[pairs_metadata[, grepl(
      pattern = "structure_taxonomy_npclassifier_",
      x = colnames(pairs_metadata)
    )]],
    sep = "|",
    direction = "long"
  ) |>
  dplyr::filter(
    !is.na(structure_taxonomy_npclassifier_01pathway) &
      !is.na(structure_taxonomy_npclassifier_02superclass) &
      !is.na(structure_taxonomy_npclassifier_03class)
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::filter(organism_taxonomy_06family %in% families_matched_restricted$unique_name) |>
  dplyr::filter(!is.na(structure_taxonomy_npclassifier_03class)) |>
  dplyr::distinct(organism_taxonomy_06family,
    structure_inchikey,
    .keep_all = TRUE
  )

specific_classes_o <- specific_classes |>
  dplyr::group_by(organism_taxonomy_06family) |>
  dplyr::distinct(organism_taxonomy_08genus, .keep_all = TRUE) |>
  dplyr::count(name = "o") %>%
  dplyr::ungroup()

tr_restricted$tip.label <-
  gsub(
    pattern = "_.*",
    replacement = "",
    x = tr_restricted$tip.label
  )

taxonomy <-
  dplyr::left_join(
    families_restricted,
    pairs_metadata
  ) |>
  dplyr::distinct(
    Domain = organism_taxonomy_01domain,
    Kingdom = organism_taxonomy_02kingdom,
    Phylum = organism_taxonomy_03phylum,
    Class = organism_taxonomy_04class,
    Order = organism_taxonomy_05order,
    Family = organism_taxonomy_06family
  )

info <- taxonomy |>
  dplyr::select(
    id = Family,
    dplyr::everything()
  ) |>
  dplyr::mutate(Kingdom = forcats::fct_reorder(Kingdom, !is.na(Domain))) |>
  dplyr::mutate(Phylum = forcats::fct_reorder(Phylum, !is.na(Kingdom))) |>
  dplyr::left_join(specific_classes_o, by = c("id" = "organism_taxonomy_06family"))

specific_classes_adapted <- specific_classes |>
  dplyr::distinct(
    structure = structure_wikidata,
    organism = organism_taxonomy_06family,
    best_candidate_1 = structure_taxonomy_npclassifier_01pathway,
    best_candidate_2 = structure_taxonomy_npclassifier_02superclass,
    best_candidate_3 = structure_taxonomy_npclassifier_03class
  ) |>
  dplyr::mutate(sample = organism, species = organism)

myHierch <-
  prepare_hierarchy(dataframe = specific_classes_adapted, type = "literature")

myPreplot <- prepare_plot(dataframe = myHierch)

ott_ready <- ott_in_tree[!duplicated(ott_in_tree) & names(ott_in_tree) %in% info$id]

tr_ott <- rotl::tol_induced_subtree(ott_ids = ott_ready)

tr_ott$tip.label <-
  gsub(
    pattern = "_.*",
    replacement = "",
    x = tr_ott$tip.label
  )

tree_ott <- ggtree::ggtree(tr = tr_ott)

p_1 <- tree_ott %<+%
  info +
  ggtree::geom_tiplab(
    ggplot2::aes(color = Kingdom),
    align = TRUE,
    size = ggplot2::rel(5),
    offset = ggplot2::rel(1)
  ) +
  ggplot2::scale_color_manual(
    values = strsplit(x = paired, split = " "),
    na.value = "grey"
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot,
    geom = geom_col,
    mapping =
      ggplot2::aes(
        y = sample,
        x = relative,
        fill = ids
      ),
    offset = ggplot2::rel(0.15),
    pwidth = ggplot2::rel(1.2),
    orientation = "y",
    stat = "identity",
  ) +
  ggplot2::scale_fill_discrete(
    name = "Chemical pathway",
    direction = "vertical",
    guide = ggplot2::guide_legend(order = 3)
  ) +
  ggplot2::scale_fill_manual(
    values = strsplit(
      x = levels(myPreplot$color) |>
        as.character(),
      split = " "
    ),
    na.value = "grey",
    guide = ggplot2::guide_legend(reverse = TRUE)
  ) +
  ggtree::geom_tippoint(mapping = ggplot2::aes(color = Kingdom, size = o)) +
  ggnewscale::new_scale_fill() +
  ggplot2::scale_fill_discrete(
    name = "Biological kingdom",
    guide = ggplot2::guide_legend(
      order = 1,
      ncol = 2
    )
  ) +
  ggplot2::scale_size_continuous(
    name = "Genera in biological family",
    guide = ggplot2::guide_legend(
      order = 2,
      direction = "horizontal"
    )
  ) +
  ggplot2::theme(
    legend.position = c(0.14, 0.93),
    legend.background = ggplot2::element_rect(fill = NA),
    legend.title = ggplot2::element_text(size = ggplot2::rel(3)),
    legend.text = ggplot2::element_text(size = ggplot2::rel(2)),
  )

p_2 <- tree_ott %<+%
  info +
  ggtree::geom_tiplab(
    ggplot2::aes(color = Kingdom),
    align = TRUE,
    size = ggplot2::rel(5),
    offset = ggplot2::rel(1)
  ) +
  ggplot2::scale_color_manual(
    values = strsplit(x = paired, split = " "),
    na.value = "grey"
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Terpenoids"),
    geom = geom_col,
    mapping =
      ggplot2::aes(
        y = sample,
        x = log(values),
        fill = ids
      ),
    offset = ggplot2::rel(0.15),
    pwidth = ggplot2::rel(0.19),
    orientation = "y",
    stat = "identity",
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Shikimates and Phenylpropanoids"),
    geom = geom_col,
    mapping =
      ggplot2::aes(
        y = sample,
        x = log(values),
        fill = ids
      ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y",
    stat = "identity",
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Alkaloids"),
    geom = geom_col,
    mapping =
      ggplot2::aes(
        y = sample,
        x = log(values),
        fill = ids
      ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y",
    stat = "identity",
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Polyketides"),
    geom = geom_col,
    mapping =
      ggplot2::aes(
        y = sample,
        x = log(values),
        fill = ids
      ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y",
    stat = "identity",
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Fatty acids"),
    geom = geom_col,
    mapping =
      ggplot2::aes(
        y = sample,
        x = log(values),
        fill = ids
      ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y",
    stat = "identity",
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Amino acids and Peptides"),
    geom = geom_col,
    mapping =
      ggplot2::aes(
        y = sample,
        x = log(values),
        fill = ids
      ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y",
    stat = "identity",
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Carbohydrates"),
    geom = geom_col,
    mapping =
      ggplot2::aes(
        y = sample,
        x = log(values),
        fill = ids
      ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.09),
    orientation = "y",
    stat = "identity",
  ) +
  ggplot2::scale_fill_discrete(
    name = "Chemical pathway",
    direction = "vertical",
    guide = ggplot2::guide_legend(order = 3)
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      strsplit(rev(unique(droplevels(myPreplot[myPreplot$parents == "Alkaloids", ]$color))) |> as.character(), split = " "),
      strsplit(rev(unique(droplevels(myPreplot[myPreplot$parents == "Amino acids and Peptides", ]$color))) |> as.character(), split = " "),
      strsplit(rev(unique(droplevels(myPreplot[myPreplot$parents == "Carbohydrates", ]$color))) |> as.character(), split = " "),
      strsplit(rev(unique(droplevels(myPreplot[myPreplot$parents == "Fatty acids", ]$color))) |> as.character(), split = " "),
      strsplit(rev(unique(droplevels(myPreplot[myPreplot$parents == "Polyketides", ]$color))) |> as.character(), split = " "),
      strsplit(rev(unique(droplevels(myPreplot[myPreplot$parents == "Shikimates and Phenylpropanoids", ]$color))) |> as.character(), split = " "),
      strsplit(rev(unique(droplevels(myPreplot[myPreplot$parents == "Terpenoids", ]$color))) |> as.character(), split = " ")
    ),
    na.value = "grey"
  ) +
  ggtree::geom_tippoint(mapping = ggplot2::aes(color = Kingdom, size = o)) +
  ggnewscale::new_scale_fill() +
  ggplot2::scale_fill_discrete(
    name = "Biological kingdom",
    guide = ggplot2::guide_legend(
      order = 1,
      ncol = 2
    )
  ) +
  ggplot2::scale_size_continuous(
    name = "Genera in biological family",
    guide = ggplot2::guide_legend(
      order = 2,
      direction = "horizontal"
    )
  ) +
  ggplot2::theme(
    legend.position = c(0.14, 0.93),
    legend.background = ggplot2::element_rect(fill = NA),
    legend.title = ggplot2::element_text(size = ggplot2::rel(3)),
    legend.text = ggplot2::element_text(size = ggplot2::rel(2)),
  )

ggplot2::ggsave(
  filename = file.path("~/Downloads/newTree_1.pdf"),
  plot = p_1,
  width = 75,
  height = 75,
  units = "in",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path("~/Downloads/newTree_2.pdf"),
  plot = p_2,
  width = 75,
  height = 75,
  units = "in",
  limitsize = FALSE
)

# test <- p_1 + p_2
#
# ggplot2::ggsave(
#   filename = file.path("~/Downloads/test.pdf"),
#   plot = test,
#   width = 150,
#   height = 75,
#   units = "in",
#   limitsize = FALSE
# )

end <- Sys.time()

log_debug("Script finished in", format(end - start))
