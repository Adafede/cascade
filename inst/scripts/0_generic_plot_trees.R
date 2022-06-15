source(file = "R/log_debug.R")
log_debug("This script plots an alternative magic tree.")

start <- Sys.time()

#' Packages
packages_cran <-
  c(
    "data.table",
    "devtools",
    "dplyr",
    "forcats",
    "ggplot2",
    "ggnewscale",
    "readr",
    "rotl",
    "splitstackshape",
    "tidyr",
    "yaml"
  )
packages_bioconductor <- c("ggtree", "ggtreeExtra", "ggstar")
packages_github <- c("KarstensLab/microshades")

source(file = "R/check_and_load_packages.R")
source(file = "R/check_export_dir.R")
source(file = "R/colors.R")
source(file = "R/load_lotus.R")
source(file = "R/make_2D.R")
source(file = "R/make_chromatographiable.R")
source(file = "R/parse_yaml_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/prepare_plot.R")

check_and_load_packages()

devtools::source_url(
  "https://raw.githubusercontent.com/taxonomicallyinformedannotation/tima-r/main/R/get_lotus.R"
)

paths <- parse_yaml_paths()
params <- parse_yaml_params()

export_name <-
  ifelse(
    test = !is.null(params$organisms$taxon),
    yes = params$organisms$taxon,
    no = "full"
  )

load_lotus()

lotus <-
  readr::read_delim(file = paths$inst$extdata$source$libraries$lotus) |>
  data.table::data.table()

if (params$structures$dimensionality == 2) {
  lotus <- lotus |>
    make_2D() |>
    data.table::data.table()
}

if (params$structures$c18 == TRUE) {
  lotus <- lotus |>
    make_chromatographiable() |>
    data.table::data.table()
}

if (!is.na(params$organisms$taxon)) {
  taxon_prerestricted <- lotus |>
    dplyr::filter(!!as.name(params$organisms$level) == params$organisms$taxon)
} else {
  taxon_prerestricted <- lotus
}

taxon_restricted <- taxon_prerestricted |>
  dplyr::filter(!is.na(!!as.name(params$organisms$group))) |>
  dplyr::distinct(structure_inchikey,
    !!as.name(params$organisms$group),
    .keep_all = TRUE
  ) |>
  dplyr::group_by(!!as.name(params$organisms$group)) |>
  dplyr::add_count() |>
  dplyr::ungroup() |>
  dplyr::filter(n >= params$structures$min) |>
  dplyr::distinct(!!as.name(params$organisms$group))

taxon_matched_restricted <-
  rotl::tnrs_match_names(
    names = taxon_restricted[[1]],
    do_approximate_matching = FALSE
  )

ott_in_tree <-
  rotl::ott_id(taxon_matched_restricted)[rotl::is_in_tree(rotl::ott_id(taxon_matched_restricted))]

taxon_restricted <- taxon_restricted |>
  dplyr::filter(!!as.name(params$organisms$group) %in% names(ott_in_tree))

taxon_matched_restricted <-
  rotl::tnrs_match_names(
    names = taxon_restricted[[1]],
    do_approximate_matching = FALSE
  )

taxon_matched_restricted <- taxon_matched_restricted |>
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

specific_classes <- taxon_prerestricted |>
  splitstackshape::cSplit(
    splitCols = colnames(lotus)[lotus[, grepl(
      pattern = "structure_taxonomy_npclassifier_",
      x = colnames(lotus)
    )]],
    sep = "|",
    direction = "long"
  ) |>
  dplyr::mutate(
    structure_taxonomy_npclassifier_01pathway = ifelse(
      test = !is.na(structure_taxonomy_npclassifier_01pathway),
      yes = structure_taxonomy_npclassifier_01pathway,
      no = "notClassified"
    )
  ) |>
  dplyr::mutate(
    structure_taxonomy_npclassifier_02superclass = ifelse(
      test = !is.na(structure_taxonomy_npclassifier_02superclass),
      yes = structure_taxonomy_npclassifier_02superclass,
      no = paste(structure_taxonomy_npclassifier_01pathway, "notClassified")
    )
  ) |>
  dplyr::mutate(
    structure_taxonomy_npclassifier_03class = ifelse(
      test = !is.na(structure_taxonomy_npclassifier_03class),
      yes = structure_taxonomy_npclassifier_03class,
      no = paste(
        structure_taxonomy_npclassifier_02superclass,
        "notClassified"
      )
    )
  ) |>
  dplyr::mutate_all(as.character) |>
  dplyr::filter(!!as.name(params$organisms$group) %in% taxon_matched_restricted$unique_name) |>
  dplyr::filter(!is.na(structure_taxonomy_npclassifier_03class)) |>
  dplyr::distinct(!!as.name(params$organisms$group),
    structure_inchikey,
    .keep_all = TRUE
  )

specific_classes_o <- specific_classes |>
  dplyr::group_by(!!as.name(params$organisms$group)) |>
  dplyr::distinct(!!as.name(params$organisms$subgroup), .keep_all = TRUE) |>
  dplyr::count(name = "o") |>
  dplyr::ungroup()

tr_restricted$tip.label <-
  gsub(
    pattern = "_.*",
    replacement = "",
    x = tr_restricted$tip.label
  )

taxonomy <-
  dplyr::left_join(
    taxon_restricted,
    lotus
  ) |>
  dplyr::distinct(
    Domain = organism_taxonomy_01domain,
    Kingdom = organism_taxonomy_02kingdom,
    Phylum = organism_taxonomy_03phylum,
    Class = organism_taxonomy_04class,
    Order = organism_taxonomy_05order,
    Family = organism_taxonomy_06family,
    Tribe = organism_taxonomy_07tribe,
    Genus = organism_taxonomy_08genus,
    Species = organism_taxonomy_09species,
    Varietas = organism_taxonomy_10varietas
  )

info <- taxonomy |>
  dplyr::select(
    id = names(taxonomy)[grepl(
      pattern = gsub(
        pattern = "organism_taxonomy_[0-9]{2}",
        replacement = "",
        x = params$organisms$group
      ),
      x = names(taxonomy),
      ignore.case = TRUE
    )],
    dplyr::everything()
  ) |>
  dplyr::mutate(Kingdom = forcats::fct_reorder(Kingdom, !is.na(Domain))) |>
  dplyr::mutate(Phylum = forcats::fct_reorder(Phylum, !is.na(Kingdom))) |>
  dplyr::left_join(specific_classes_o, by = c("id" = params$organisms$group)) |>
  dplyr::mutate(id = gsub(
    pattern = " ",
    replacement = "_",
    x = id
  ))

specific_classes_adapted <- specific_classes |>
  dplyr::distinct(
    structure = structure_wikidata,
    organism = !!as.name(params$organisms$group),
    best_candidate_1 = structure_taxonomy_npclassifier_01pathway,
    best_candidate_2 = structure_taxonomy_npclassifier_02superclass,
    best_candidate_3 = structure_taxonomy_npclassifier_03class
  ) |>
  dplyr::mutate(sample = organism, species = organism) |>
  dplyr::mutate(
    best_candidate_1 = ifelse(
      test = best_candidate_1 == "notClassified",
      yes = "Other",
      no = best_candidate_1
    )
  )

myHierch <-
  prepare_hierarchy(dataframe = specific_classes_adapted, type = "literature")

myPreplot <- prepare_plot(dataframe = myHierch) |>
  mutate(sample = gsub(
    pattern = " ",
    replacement = "_",
    x = sample
  ))

ott_ready <-
  ott_in_tree[!duplicated(ott_in_tree) &
    gsub(
      pattern = " ",
      replacement = "_",
      x = names(ott_in_tree)
    ) %in% info$id]

tr_ott <- rotl::tol_induced_subtree(ott_ids = ott_ready)

tr_ott$tip.label <-
  gsub(
    pattern = "_ott.*",
    replacement = "",
    x = tr_ott$tip.label
  )

tree_ott <- ggtree::ggtree(tr = tr_ott)

p <- tree_ott %<+%
  info +
  ggtree::geom_tiplab(
    mapping = ggplot2::aes(color = Kingdom),
    align = TRUE,
    size = ggplot2::rel(5),
    offset = ggplot2::rel(0.2)
  ) +
  ggplot2::scale_color_manual(
    values = strsplit(
      x = c(
        nice_colors[[5]][5],
        nice_colors[[6]][5],
        nice_colors[[7]][5],
        nice_colors[[8]][5],
        nice_colors[[9]][5],
        nice_colors[[10]][5]
      ),
      split = " "
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
  ggnewscale::new_scale("size") +
  ggplot2::scale_size_continuous(
    name = "Considered subtaxa in taxon",
    guide = ggplot2::guide_legend(
      order = 2,
      direction = "horizontal"
    )
  ) +
  ggplot2::theme(
    legend.position = c(ggplot2::rel(0.15), ggplot2::rel(0.93)),
    legend.background = ggplot2::element_rect(fill = NA),
    legend.title = ggplot2::element_text(size = ggplot2::rel(3)),
    legend.text = ggplot2::element_text(size = ggplot2::rel(2)),
  ) +
  ggnewscale::new_scale_fill()

p_rel <- p +
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
  ) +
  ggplot2::scale_fill_manual(
    values = strsplit(
      x = levels(myPreplot$color) |>
        as.character(),
      split = " "
    ),
    na.value = "grey",
    guide = ggplot2::guide_legend(reverse = TRUE)
  )

p_abs <- p +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Terpenoids"),
    geom = geom_col,
    mapping = ggplot2::aes(
      y = sample,
      x = log(values),
      fill = ids
    ),
    offset = ggplot2::rel(0.15),
    pwidth = ggplot2::rel(0.19),
    orientation = "y"
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Shikimates and Phenylpropanoids"),
    geom = geom_col,
    mapping = ggplot2::aes(
      y = sample,
      x = log(values),
      fill = ids
    ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y"
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Alkaloids"),
    geom = geom_col,
    mapping = ggplot2::aes(
      y = sample,
      x = log(values),
      fill = ids
    ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y"
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
    orientation = "y"
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Fatty acids"),
    geom = geom_col,
    mapping = ggplot2::aes(
      y = sample,
      x = log(values),
      fill = ids
    ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y"
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Amino acids and Peptides"),
    geom = geom_col,
    mapping = ggplot2::aes(
      y = sample,
      x = log(values),
      fill = ids
    ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.19),
    orientation = "y"
  ) +
  ggtreeExtra::geom_fruit(
    data = myPreplot |> dplyr::filter(parents == "Carbohydrates"),
    geom = geom_col,
    mapping = ggplot2::aes(
      y = sample,
      x = log(values),
      fill = ids
    ),
    offset = ggplot2::rel(0.0),
    pwidth = ggplot2::rel(0.09),
    orientation = "y"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      strsplit(rev(unique(
        droplevels(myPreplot[myPreplot$parents == "Alkaloids", ]$color)
      )) |> as.character(), split = " "),
      strsplit(rev(unique(
        droplevels(myPreplot[myPreplot$parents == "Amino acids and Peptides", ]$color)
      )) |> as.character(), split = " "),
      strsplit(rev(unique(
        droplevels(myPreplot[myPreplot$parents == "Carbohydrates", ]$color)
      )) |> as.character(), split = " "),
      strsplit(rev(unique(
        droplevels(myPreplot[myPreplot$parents == "Fatty acids", ]$color)
      )) |> as.character(), split = " "),
      strsplit(rev(unique(
        droplevels(myPreplot[myPreplot$parents == "Polyketides", ]$color)
      )) |> as.character(), split = " "),
      strsplit(rev(unique(
        droplevels(myPreplot[myPreplot$parents == "Shikimates and Phenylpropanoids", ]$color)
      )) |> as.character(), split = " "),
      strsplit(rev(unique(
        droplevels(myPreplot[myPreplot$parents == "Terpenoids", ]$color)
      )) |> as.character(), split = " ")
    ),
    na.value = "grey"
  )

lapply(X = paths$data$trees$path, FUN = check_export_dir)

ggplot2::ggsave(
  filename = file.path(
    paths$data$trees$path,
    paste("tree", export_name, "relative.pdf", sep = "_")
  ),
  plot = p_rel,
  width = 75,
  height = 70 + nrow(taxon_restricted) / 10,
  units = "in",
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = file.path(
    paths$data$trees$path,
    paste("tree", export_name, "absolute.pdf", sep = "_")
  ),
  plot = p_abs,
  width = 75,
  height = 70 + nrow(taxon_restricted) / 10,
  units = "in",
  limitsize = FALSE
)

end <- Sys.time()

log_debug("Script finished in", format(end - start))
