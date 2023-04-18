#' paths
path_peaks_pos <-
  "data/interim/peaks/191109_AR_10043_enriched_UHR_Pos_featuresInformed_cad.tsv.gz"
path_peaks_neg <-
  "data/interim/peaks/191109_AR_10043_enriched_UHR_Neg_featuresInformed_cad.tsv.gz"
path_annotations_final <- "data/paper/final_table.tsv"
path_lotus <- "data/source/libraries/lotus.csv.gz"

source(file = "R/get_params.R")
source(file = "R/parse_cli_params.R")
source(file = "R/parse_yaml_paths.R")
source(file = "R/parse_yaml_params.R")
step <- "processing"
paths <- parse_yaml_paths()
params <- ""
params <- get_params(step = step)

#' import
tables_peaks_pos <- readr::read_delim(
  file = path_peaks_pos,
  col_select = c("peak_rt_apex", "peak_area")
) |>
  dplyr::distinct()
table_peaks_neg <- readr::read_delim(
  file = path_peaks_neg,
  col_select = c("peak_rt_apex", "peak_area")
) |>
  dplyr::distinct()
table_annotations_initial <-
  readr::read_delim(file = path_annotations_final)
lotus <- readr::read_delim(
  file = path_lotus,
  col_select = c(
    "structure_inchikey",
    "organism_taxonomy_01domain",
    "organism_taxonomy_02kingdom",
    "organism_taxonomy_03phylum",
    "organism_taxonomy_04class",
    "organism_taxonomy_05order",
    "organism_taxonomy_06family",
    "organism_taxonomy_07tribe",
    "organism_taxonomy_08genus",
    "organism_taxonomy_09species",
    "organism_taxonomy_10varietas",
    "reference_doi"
  )
) |>
  dplyr::distinct()

#' process
table_peaks <- tables_peaks_pos |>
  dplyr::bind_rows(table_peaks_neg) |>
  dplyr::mutate(
    rt_min = peak_rt_apex - params$chromato$peak$tolerance$rt / 2,
    rt_max = peak_rt_apex + params$chromato$peak$tolerance$rt / 2
  ) |>
  data.table::data.table()

data.table::setkey(table_peaks, rt_min, rt_max)

table_processed <- table_annotations_initial |>
  dplyr::select(
    peak_rt_apex,
    peak_area,
    feature_id,
    feature_rt = rt,
    feature_mz = mz,
    smiles_2D,
    best_candidate_1,
    best_candidate_2,
    best_candidate_3,
    score_final,
    best_candidate_organism,
    inchikey_2D
  ) |>
  dplyr::mutate(
    rt_min = peak_rt_apex,
    rt_max = peak_rt_apex
  ) |>
  data.table::data.table()

data.table::setkey(table_processed, rt_min, rt_max)

log_debug(x = "joining within given rt tolerance")
table_medium <-
  data.table::foverlaps(table_peaks, table_processed) |>
  dplyr::arrange(dplyr::desc(score_final)) |>
  dplyr::distinct(i.peak_rt_apex, smiles_2D, .keep_all = TRUE) |>
  dplyr::mutate(rt_min = ifelse(
    test = !is.na(rt_min),
    yes = rt_min,
    no = round(x = i.peak_rt_apex, digits = 1)
  )) |>
  dplyr::group_by(rt_min) |>
  dplyr::mutate(
    newrt = mean(i.peak_rt_apex),
    newarea = mean(i.peak_area)
  ) |>
  dplyr::ungroup() |>
  dplyr::distinct(newrt, newarea, .keep_all = TRUE) |>
  dplyr::arrange(newrt) |>
  dplyr::mutate(diff = newrt - dplyr::lag(newrt)) |> # eliminate rounding duplicates
  dplyr::filter(!is.na(diff) &
    diff > params$chromato$peak$tolerance$rt / 2) |>
  dplyr::mutate(peak_id = dplyr::row_number()) |>
  dplyr::select(
    peak_id,
    peak_rt = newrt,
    peak_area = newarea,
    feature_id,
    feature_rt,
    feature_mz,
    smiles_2D,
    best_candidate_1,
    best_candidate_2,
    best_candidate_3,
    score_final,
    best_candidate_organism,
    inchikey_2D
  ) |>
  dplyr::mutate(
    peak_rt = round(x = peak_rt, digits = 2),
    peak_area = round(x = 100 * peak_area / sum(peak_area), digits = 2),
    feature_rt = round(x = feature_rt, digits = 2),
    feature_mz = round(x = feature_rt, digits = 4),
    score_final = round(x = score_final, digits = 2)
  )

lotus_long <- lotus |>
  dplyr::mutate(inchikey_2D = gsub(
    pattern = "-.*",
    replacement = "",
    x = structure_inchikey
  )) |>
  tidyr::pivot_longer(cols = 2:11, names_prefix = "organism_taxonomy_") |>
  dplyr::distinct(inchikey_2D,
    reference = reference_doi,
    best_candidate_organism = value
  )

table_final <- table_medium |>
  dplyr::left_join(lotus_long) |>
  dplyr::group_by(peak_id, inchikey_2D) |>
  dplyr::mutate(reference = paste(reference, collapse = "; ")) |>
  dplyr::mutate(Structure = RCurl::curlEscape(smiles_2D)) |>
  dplyr::mutate(
    reference = ifelse(
      test = reference == "NA",
      yes = NA,
      no = reference
    ),
    Structure = ifelse(
      test = Structure == "NA",
      yes = "",
      no = Structure
    ),
    score_final = ifelse(
      test = score_final > 1,
      yes = as.character(">1.00"),
      no = as.character(score_final)
    ),
    best_candidate_organism = gsub(
      pattern = "Swertia chirayita",
      replacement = "S. chirayita",
      x = best_candidate_organism
    ),
    best_candidate_1 = ifelse(
      test = best_candidate_1 == "Shikimates and Phenylpropanoids",
      yes = "Shik \\& PAL",
      no = best_candidate_1
    ),
    best_candidate_2 = ifelse(
      test = best_candidate_2 == "Fatty Acids and Conjugates",
      yes = "FA \\& conj.",
      no = best_candidate_2
    ),
    best_candidate_2 = ifelse(
      test = best_candidate_2 == "Polycyclic aromatic polyketides",
      yes = "Polyc. aro. polyk.",
      no = best_candidate_2
    ),
    best_candidate_2 = ifelse(
      test = best_candidate_2 == "Phenylethanoids (C6-C2)",
      yes = "$\\phi$ eth. (C6-C2)",
      no = best_candidate_2
    ),
    best_candidate_3 = best_candidate_3 |>
      gsub(pattern = "sesquiterpenoids", replacement = "sesquiter.") |>
      gsub(pattern = "triterpenoids", replacement = "triter.") |>
      gsub(pattern = "diterpenoids", replacement = "diter.") |>
      gsub(pattern = "monoterpenoids", replacement = "monoter.") |>
      gsub(pattern = " and ", replacement = " \\& ", fixed = TRUE) |>
      gsub(pattern = "Phyllocladane", replacement = "Phyllocl.") |>
      gsub(pattern = "Anthraquinones", replacement = "Anthraq.") |>
      gsub(pattern = "Taraxastane", replacement = "Taraxa.") |>
      gsub(pattern = "Kaurane", replacement = "Kaur.") |>
      gsub(pattern = "Ursane", replacement = "Urs.") |>
      gsub(pattern = "Moretane", replacement = "Moret.")
  ) |>
  dplyr::ungroup() |>
  dplyr::distinct() |>
  dplyr::select(
    `Peak ID` = peak_id,
    `Peak RT [min]` = peak_rt,
    `Peak Area [%]` = peak_area,
    `Feature ID` = feature_id,
    `Feature RT [min]` = feature_rt,
    `Feature m/z` = feature_mz,
    Structure,
    `InChIKey 2D` = inchikey_2D,
    `Annotation Score` = score_final,
    `Chemical Pathway` = best_candidate_1,
    `Chemical Superclass` = best_candidate_2,
    `Chemical Class` = best_candidate_3,
    `Closest Organism with Structure Reported` = best_candidate_organism,
    `Reference(s)` = reference
  )

test <- table_final |>
  dplyr::mutate_if(is.character, ~ replace(., is.na(.), "Unknown")) |>
  dplyr::distinct(
    `Peak ID`,
    `Peak Area [%]`,
    `Chemical Class`
  ) |>
  gt::gt(rowname_col = "row", groupname_col = "Chemical Class") |>
  gt::summary_rows(
    groups = TRUE,
    columns = `Peak Area [%]`,
    fns = list(Total = "sum")
  ) |>
  gt::sub_missing(
    columns = everything(),
    rows = everything(),
    missing_text = ""
  )
test

source("r/molinfo.R")

prettyTable <- table_final |>
  gt::gt(auto_align = TRUE) |>
  # cols_width(
  #   everything() ~ px(200)
  # ) |>
  gt::text_transform(
    locations = gt::cells_body(columns = Structure),
    fn = molinfo
  ) |>
  gt::sub_missing(
    columns = everything(),
    rows = everything(),
    missing_text = ""
  )

prettyTable

table_final_2 <- table_final |>
  dplyr::rowwise() |>
  dplyr::mutate(Structure = ifelse(
    test = !is.na(`InChIKey 2D`),
    yes = paste0("\\includegraphics[width=0.33333in,height=0.075in]{images/", `InChIKey 2D`, ".pdf}"),
    no = NA
  ))
readr::write_tsv(x = table_final_2, file = "data/paper/prettyTable.tsv", na = "")
gt::gtsave(
  data = prettyTable,
  filename = "data/paper/prettyTable.html"
)
gt::gtsave(
  data = test,
  filename = "data/paper/prettyTable_2.html"
)
