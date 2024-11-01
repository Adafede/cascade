start <- Sys.time()

source(file = "R/add_peak_metadata.R")
source(file = "R/keep_best_candidates.R")
source(file = "R/make_confident.R")
source(file = "R/make_other.R")
source(file = "R/no_other.R")
source(file = "R/plot_peaks_statistics.R")
source(file = "R/plot_results.R")
source(file = "R/prepare_comparison.R")
source(file = "R/prepare_hierarchy.R")
source(file = "R/preprocess_chromatograms.R")
source(file = "R/treemaps_progress.R")
source(file = "R/cascade-package.R")

message(
  "This program performs",
  "TODO"
)
message("Authors: \n", "AR")
message("Contributors: \n", "...")

#' Specific paths
ANNOTATIONS <- "~/.tima/data/processed/241026_103144_extract/extract_results.tsv"
BPI <- FALSE
CAD <- TRUE
PDA <- FALSE
CAD_SHIFT <- 0.05
PDA_SHIFT <- 0.1
TIME_MIN <- 0.7
TIME_MAX <- 35.2
CONFIDENCE_SCORE_MIN <- 0.4
PEAK_SIMILARITY_PREFILTER <- 0.6
PEAK_SIMILARITY <- 0.8
THESIS <- FALSE
IMPORT_FILE_CAD <- "data/interim/peaks/210619_AR_06_V_03_2_01_featuresInformed_cad.tsv.gz"
IMPORT_FILE_CAD_2 <- "data/interim/peaks/210619_AR_06_V_03_2_01_featuresInformed_cad.tsv.gz"
EXPORT_DIR <- "data/interim/peaks"
FILE_POSITIVE <- "data/source/mzml/210619_AR_06_V_03_2_01.mzML"
names <- FILE_POSITIVE |>
  gsub(pattern = ".*/", replacement = "") |>
  gsub(pattern = "[0-9]{8}_AR_[0-9]{2}_", replacement = "") |>
  gsub(
    pattern = ".mzML",
    replacement = "",
    fixed = TRUE
  )

message("listing files")
message("loading raw files (can take long if loading multiple files)")
dda_data <- MSnbase::readMSData(
  files = FILE_POSITIVE,
  mode = "onDisk",
  msLevel. = 1
)

if (THESIS == TRUE) {
  dda_data_neg <-
    MSnbase::readMSData(
      files = files |>
        gsub(pattern = "_Pos", replacement = "_Neg") |>
        gsub(pattern = "191113", replacement = "191107"),
      mode = "onDisk",
      msLevel. = 1
    )
}

message("opening raw files objects and extracting chromatograms")
chromatograms_all <- lapply(FILE_POSITIVE, mzR::openMSfile) |>
  lapply(mzR::chromatograms) |>
  purrr::flatten()

if (THESIS == TRUE) {
  chromatograms_all_neg <-
    lapply(
      files |> gsub(pattern = "_Pos", replacement = "_Neg") |>
        gsub(pattern = "191113", replacement = "191107"),
      mzR::openMSfile
    ) |>
    lapply(mzR::chromatograms) |>
    purrr::flatten()
}

chromatograms_list_bpi <- preprocess_chromatograms(
  detector = "bpi",
  list = chromatograms_all[c(TRUE, FALSE, FALSE)],
  signal_name = "BasePeak_0",
  shift = 0
)
chromatograms_list_cad <- preprocess_chromatograms()
chromatograms_list_pda <-
  preprocess_chromatograms(
    detector = "pda",
    list = chromatograms_all[c(FALSE, TRUE, FALSE)],
    signal_name = "PDA.1_TotalAbsorbance_0",
    shift = PDA_SHIFT
  )

if (THESIS == TRUE) {
  chromatograms_list_bpi_neg <- preprocess_chromatograms(
    detector = "bpi",
    list = chromatograms_all_neg[c(TRUE, FALSE, FALSE)],
    signal_name = "BasePeak_0",
    shift = 0
  )
  chromatograms_list_cad_neg <-
    preprocess_chromatograms(list = chromatograms_all_neg[c(FALSE, FALSE, TRUE)])
  chromatograms_list_pda_neg <-
    preprocess_chromatograms(
      detector = "pda",
      list = chromatograms_all_neg[c(FALSE, TRUE, FALSE)],
      signal_name = "PDA.1_TotalAbsorbance_0",
      shift = PDA_SHIFT
    )
}

message("loading annotations")
annotations <-
  ANNOTATIONS |>
  lapply(
    FUN = function(x) {
      readr::read_delim(file = x) |>
        dplyr::mutate(dplyr::across(
          dplyr::contains("candidate_structure_tax"),
          .fns = function(x) {
            gsub(
              pattern = "\\$",
              replacement = "or",
              x = x
            )
          }
        )) |>
        dplyr::mutate(dplyr::across(
          dplyr::contains("candidate_structure_tax"),
          .fns = function(x) {
            gsub(
              pattern = "§",
              replacement = "$",
              x = x
            )
          }
        )) |>
        dplyr::mutate(mode = ifelse(
          test = grepl(
            pattern = "_pos.",
            x = x,
            ignore.case = TRUE
          ),
          yes = "pos",
          no = "neg"
        ))
    }
  ) |>
  dplyr::bind_rows()

message("keeping best annotations only")
best_candidates <- annotations |>
  keep_best_candidates()

message("adding metadata dirtily for now")
candidates_metadata <- best_candidates |>
  dplyr::mutate(species = "Swertia chirayita") |>
  dplyr::mutate(feature_id = as.numeric(feature_id))

message("keeping only candidates above desired threshold")
candidates_confident <- candidates_metadata |>
  make_confident(score = CONFIDENCE_SCORE_MIN)

myDirtyListF <- function(list, mode = "pos") {
  purrr::map(.x = list, .f = ~ dplyr::filter(.x, grepl(
    pattern = paste0("_", !!as.character(mode)),
    x = sample,
    ignore.case = TRUE
  )))
}

if (BPI) {
  detector <- "bpi"
}
if (BPI) {
  compared_peaks_list_bpi <- prepare_comparison(detector = "bpi")
  compared_peaks_list_bpi_pos <- compared_peaks_list_bpi |>
    myDirtyListF()
  compared_peaks_list_bpi_neg <- compared_peaks_list_bpi |>
    myDirtyListF(mode = "neg")
  plots_1_bpi_pos <- plot_results_1(detector = "bpi")
  plots_1_bpi_neg <- plot_results_1(detector = "bpi", mode = "neg")
  plots_2_bpi <- plot_results_2(detector = "bpi")
}
if (CAD) {
  detector <- "cad"
}
# TODO investigate
mode <- "pos"

if (CAD) {
  compared_peaks_list_cad <- prepare_comparison()
  ## TODO FIX THIS
  compared_peaks_list_cad_pos <- compared_peaks_list_cad
  compared_peaks_list_cad_neg <- compared_peaks_list_cad |>
    myDirtyListF(mode = "neg")
  plots_1_cad_pos <- plot_results_1()
  plots_1_cad_neg <- plot_results_1(mode = "neg")
  plots_2_cad <- plot_results_2()
}

if (PDA) {
  detector <- "pda"
}

if (PDA) {
  compared_peaks_list_pda <- prepare_comparison(detector = "pda")
  compared_peaks_list_pda_pos <- compared_peaks_list_pda |>
    myDirtyListF()
  compared_peaks_list_pda_neg <- compared_peaks_list_pda |>
    myDirtyListF(mode = "neg")
  plots_1_pda_pos <- plot_results_1(detector = "pda")
  plots_1_pda_neg <- plot_results_1(detector = "pda", mode = "neg")
  plots_2_pda <- plot_results_2(detector = "pda")
}

fig_taxo <-
  ggpubr::ggarrange(
    plots_1_cad_pos$histograms_taxo_maj,
    ncol = 1,
    nrow = 2,
    plots_1_cad_neg$histograms_taxo_maj,
    align = "hv",
    common.legend = TRUE,
    legend = "bottom",
    labels = "AUTO"
  )
fig_taxo

fig_minmaj <-
  ggpubr::ggarrange(
    plots_1_cad_pos$histograms_unique_conf_maj,
    ncol = 2,
    nrow = 1,
    plots_1_cad_pos$histograms_unique_conf_min,
    align = "hv",
    common.legend = TRUE,
    legend = "bottom"
  )
fig_minmaj

# test_taxo_3 <-
#   ggpubr::ggarrange(
#     ncol = 2,
#     nrow = 2,
#     plots_1_cad_pos$histograms_unique_conf_maj,
#     plots_1_cad_neg$histograms_unique_conf_maj,
#     plots_1_cad_pos$histograms_unique_conf_min,
#     plots_1_cad_neg$histograms_unique_conf_min,
#     common.legend = TRUE,
#     legend = "bottom",
#     labels = "AUTO"
#   )
# test_taxo_3

hierarchies <- list()
hierarchies$peaks_maj <-
  prepare_hierarchy(
    dataframe = compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
      # dplyr::filter(id == "bitter/191113_AR_10043_Pos") |>
      # dplyr::filter(id == "bitter/TODO") |>
      dplyr::filter(id == "UHR/191109_AR_10043_enriched_UHR_Pos") |>
      dplyr::mutate(id = "peaks_maj") |>
      dplyr::mutate(sample = id, species = id) |>
      dplyr::distinct() |>
      make_other() |>
      no_other(),
    type = "analysis",
    detector = "cad"
  )
hierarchies$peaks_min <-
  prepare_hierarchy(
    dataframe = compared_peaks_list_cad$peaks_min_precor_taxo_cor |>
      dplyr::filter(mode == "pos") |>
      dplyr::mutate(id = "peaks_min") |>
      dplyr::mutate(sample = id, species = id) |>
      dplyr::distinct() |>
      make_other() |>
      no_other() |>
      dplyr::mutate(intensity = 1),
    # because MS intensity not OK
    type = "analysis",
    detector = "ms"
  )
hierarchies$special <- prepare_hierarchy(
  dataframe = compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
    # dplyr::filter(id == "bitter/191113_AR_10043_Pos") |>
    # dplyr::filter(id == "bitter/TODO") |>
    dplyr::filter(id == "UHR/191109_AR_10043_enriched_UHR_Pos") |>
    dplyr::filter(mode == "pos") |>
    dplyr::mutate(id = "peaks_maj") |>
    dplyr::mutate(sample = id, species = id) |>
    dplyr::select(-taxo, -taxo_2, -sum, -sum_2, -keep) |>
    rbind(
      compared_peaks_list_cad$peaks_min_precor_taxo_cor |>
        dplyr::filter(id == "UHR/191109_AR_10043_enriched_UHR_Pos") |>
        dplyr::filter(mode == "pos") |>
        dplyr::mutate(id = "peaks_min") |>
        dplyr::mutate(sample = id, species = id)
    ) |>
    dplyr::distinct(),
  type = "analysis",
  detector = "cad"
) |>
  dplyr::arrange(sample)
treemaps <-
  treemaps_progress_no_title(xs = names(hierarchies)[!grepl(
    pattern = "_grouped",
    x = names(hierarchies)
  )])
# treemaps$special

if (BPI) {
  df_meta_bpi_pos <- compared_peaks_list_bpi$peaks_all |>
    dplyr::filter(mode == "pos") |>
    dplyr::full_join(
      best_candidates |>
        dplyr::filter(mode == "pos") |>
        dplyr::mutate(feature_id = as.numeric(feature_id))
    ) |>
    dplyr::left_join(
      compared_peaks_list_bpi$peaks_maj_precor_taxo_cor |>
        dplyr::distinct(sample, peak_id, mode, feature_id, keep)
    ) |>
    add_peak_metadata()
  df_meta_bpi_neg <- compared_peaks_list_bpi$peaks_all |>
    dplyr::filter(mode == "neg") |>
    dplyr::full_join(
      best_candidates |>
        dplyr::filter(mode == "neg") |>
        dplyr::mutate(feature_id = as.numeric(feature_id))
    ) |>
    dplyr::left_join(
      compared_peaks_list_bpi$peaks_maj_precor_taxo_cor |>
        dplyr::distinct(sample, peak_id, mode, feature_id, keep)
    ) |>
    add_peak_metadata()
}

if (CAD) {
  df_meta_cad_pos <- compared_peaks_list_cad$peaks_all |>
    dplyr::filter(mode == "pos") |>
    dplyr::full_join(
      best_candidates |>
        dplyr::filter(mode == "pos") |>
        dplyr::mutate(feature_id = as.numeric(feature_id))
    ) |>
    dplyr::left_join(
      compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
        dplyr::distinct(sample, peak_id, mode, feature_id, keep)
    ) |>
    add_peak_metadata()
  df_meta_cad_neg <- compared_peaks_list_cad$peaks_all |>
    dplyr::filter(mode == "neg") |>
    dplyr::full_join(
      best_candidates |>
        dplyr::filter(mode == "neg") |>
        dplyr::mutate(feature_id = as.numeric(feature_id))
    ) |>
    dplyr::left_join(
      compared_peaks_list_cad$peaks_maj_precor_taxo_cor |>
        dplyr::distinct(sample, peak_id, mode, feature_id, keep)
    ) |>
    add_peak_metadata()
}

if (PDA) {
  df_meta_pda_pos <- compared_peaks_list_pda$peaks_all |>
    dplyr::filter(mode == "pos") |>
    dplyr::full_join(
      best_candidates |>
        dplyr::filter(mode == "pos") |>
        dplyr::mutate(feature_id = as.numeric(feature_id))
    ) |>
    dplyr::left_join(
      compared_peaks_list_pda$peaks_maj_precor_taxo_cor |>
        dplyr::distinct(sample, peak_id, mode, feature_id, keep)
    ) |>
    add_peak_metadata()
  df_meta_pda_neg <- compared_peaks_list_pda$peaks_all |>
    dplyr::filter(mode == "neg") |>
    dplyr::full_join(
      best_candidates |>
        dplyr::filter(mode == "neg") |>
        dplyr::mutate(feature_id = as.numeric(feature_id))
    ) |>
    dplyr::left_join(
      compared_peaks_list_pda$peaks_maj_precor_taxo_cor |>
        dplyr::distinct(sample, peak_id, mode, feature_id, keep)
    ) |>
    add_peak_metadata()
}

if (BPI) {
  plots_bpi_pos <- df_meta_bpi_pos |>
    plot_peaks_statistics()
  plots_bpi_neg <- df_meta_bpi_neg |>
    plot_peaks_statistics()
}

if (CAD) {
  plots_cad_pos <- df_meta_cad_pos |>
    plot_peaks_statistics()
  plots_cad_neg <- df_meta_cad_neg |>
    plot_peaks_statistics()
}

if (PDA) {
  plots_pda_pos <- df_meta_pda_pos |>
    plot_peaks_statistics()
  plots_pda_neg <- df_meta_pda_neg |>
    plot_peaks_statistics()
}

example_peak <- compared_peaks_list_cad[["peaks_all"]] |>
  dplyr::filter(!is.na(peak_id)) |>
  dplyr::mutate(newrt = round(peak_rt_apex, 1)) |>
  dplyr::select(
    peak_id,
    peak_rt_min,
    peak_rt_max,
    newrt,
    feature_id,
    feature_rt,
    feature_mz,
    comparison_score,
    mode
  ) |>
  dplyr::group_by(newrt) |>
  dplyr::add_count(sort = TRUE) |>
  dplyr::ungroup() |>
  dplyr::filter(n == unique(n)[3]) |>
  dplyr::filter(newrt == max(newrt)) # avoid equal n

list_mz_pos <- example_peak |>
  dplyr::filter(mode == "pos") |>
  dplyr::pull(feature_mz) |>
  as.list()
list_mz_neg <- example_peak |>
  dplyr::filter(mode == "neg") |>
  dplyr::pull(feature_mz) |>
  as.list()
rt_min <- example_peak |>
  dplyr::filter(peak_rt_min == min(peak_rt_min)) |>
  dplyr::distinct(peak_rt_min) |>
  dplyr::pull(peak_rt_min)
rt_max <- example_peak |>
  dplyr::filter(peak_rt_max == max(peak_rt_max)) |>
  dplyr::distinct(peak_rt_max) |>
  dplyr::pull(peak_rt_max)

aa_pos <- lapply(
  X = seq_along(list_mz_pos),
  FUN = function(x) {
    mzR::chromatogram(
      object = dda_data,
      rt = c((rt_min) * 60, (rt_max) * 60),
      mz = list_mz_pos[[x]] + c(0.01, -0.01)
    )
  }
)
aa_neg <- lapply(
  X = seq_along(list_mz_neg),
  FUN = function(x) {
    mzR::chromatogram(
      object = dda_data_neg,
      rt = c((rt_min) * 60, (rt_max) * 60),
      mz = list_mz_neg[[x]] + c(0.01, -0.01)
    )
  }
)
bb_pos <- lapply(
  X = seq_along(aa_pos),
  FUN = function(x) {
    time <- aa_pos[[x]][[1]]@rtime / 60
    intensity <-
      aa_pos[[x]][[1]]@intensity # / max(aa_pos[[x]][[1]]@intensity)
    mz <- round(aa_pos[[x]][[1]]@mz, 1)
    return(data.frame(time, intensity, mz = mean(mz)))
  }
)
bb_neg <- lapply(
  X = seq_along(aa_neg),
  FUN = function(x) {
    time <- aa_neg[[x]][[1]]@rtime / 60
    intensity <-
      aa_neg[[x]][[1]]@intensity # / max(aa_neg[[x]][[1]]@intensity)
    mz <- round(aa_neg[[x]][[1]]@mz, 1)
    return(data.frame(time, intensity, mz = mean(mz)))
  }
)

temp_df <- example_peak |>
  dplyr::mutate(mz = round(feature_mz, 1)) |>
  dplyr::distinct(feature_id, comparison_score, mz) |>
  dplyr::mutate(comparison_score = ifelse(test = comparison_score < 0,
    yes = 0,
    no = comparison_score
  )) |>
  dplyr::left_join(best_candidates |>
    dplyr::select(-mz) |>
    dplyr::mutate(feature_id = as.numeric(feature_id))) |>
  dplyr::group_by(molecular_formula) |>
  dplyr::add_count(name = "mf") |>
  dplyr::group_by(inchikey_2D) |>
  dplyr::add_count(name = "ik") |>
  dplyr::rowwise() |>
  dplyr::mutate(
    molecular_formula = ifelse(
      test = mf >= 2 &
        comparison_score >= 0.5 &
        score_biological >= 0.6,
      yes = molecular_formula,
      no = "other"
    ),
    inchikey_2D = ifelse(
      test = ik >= 2 &
        comparison_score >= 0.5 &
        score_biological >= 0.6,
      yes = inchikey_2D,
      no = "other"
    )
  ) |>
  dplyr::ungroup()

cc_pos <- dplyr::bind_rows(bb_pos) |>
  dplyr::filter(!is.na(intensity)) |>
  dplyr::filter(time >= rt_min & time <= rt_max) |>
  dplyr::mutate(intensity = intensity / max(intensity)) |>
  dplyr::mutate(mode = "pos") |>
  dplyr::left_join(temp_df) |>
  dplyr::filter(!is.na(inchikey_2D))
cc_neg <- dplyr::bind_rows(bb_neg) |>
  dplyr::filter(!is.na(intensity)) |>
  dplyr::filter(time >= rt_min & time <= rt_max) |>
  dplyr::mutate(intensity = -intensity / max(intensity)) |>
  dplyr::mutate(mode = "neg") |>
  dplyr::left_join(temp_df) |>
  dplyr::filter(!is.na(inchikey_2D))

cc_cad_pos <- chromatograms_list_cad$chromatograms_improved_long |>
  # cc_cad_pos <- chromatograms_list_cad$chromatograms_original_long |>
  dplyr::mutate(time = time + CAD_SHIFT) |>
  dplyr::filter(time >= rt_min & time <= rt_max) |>
  dplyr::mutate(
    mz = "cad signal",
    comparison_score = "cad signal",
    molecular_formula = "cad signal",
    inchikey_2D = "cad signal",
    score_biological = "cad signal"
  ) |>
  dplyr::mutate(intensity = intensity / max(intensity))

cc_cad_neg <-
  chromatograms_list_cad_neg$chromatograms_improved_long |>
  # cc_cad_neg <- chromatograms_list_cad_neg$chromatograms_original_long |>
  dplyr::mutate(time = time + CAD_SHIFT) |>
  dplyr::filter(time >= rt_min & time <= rt_max) |>
  dplyr::mutate(
    mz = "cad signal",
    comparison_score = "cad signal",
    molecular_formula = "cad signal",
    inchikey_2D = "cad signal",
    score_biological = "cad signal"
  ) |>
  dplyr::mutate(intensity = -intensity / max(intensity))

plot_comparison <- ggplot2::ggplot(
  data = NULL,
  mapping = ggplot2::aes(
    x = time,
    y = intensity,
    group = mz,
    # color = as.character(inchikey_2D)
    color = comparison_score
  )
) +
  ggplot2::geom_line(data = cc_pos) +
  ggplot2::geom_line(data = cc_neg) +
  ggplot2::geom_line(
    data = cc_cad_pos,
    color = "black",
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    data = cc_cad_neg,
    color = "black",
    linetype = "dashed"
  ) +
  viridis::scale_color_viridis(
    option = "D",
    name = "Peak Similarity Score"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    complete = FALSE,
    # legend.position = "none",
    validate = TRUE
  ) +
  xlab("Time [min]") +
  ylab("Normalized Intensity")

plot_taxo <- ggplot2::ggplot(
  data = NULL,
  mapping = ggplot2::aes(
    x = time,
    y = intensity,
    group = mz,
    # color = as.character(inchikey_2D)
    color = as.numeric(score_biological)
  )
) +
  ggplot2::geom_line(data = cc_pos) +
  ggplot2::geom_line(data = cc_neg) +
  ggplot2::geom_line(
    data = cc_cad_pos,
    color = "black",
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    data = cc_cad_neg,
    color = "black",
    linetype = "dashed"
  ) +
  viridis::scale_color_viridis(
    option = "D",
    name = "Taxonomic Distance Score"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    complete = FALSE,
    # legend.position = "none",
    validate = TRUE
  ) +
  xlab("Time [min]") +
  ylab("Normalized Intensity")

plot_mf <- ggplot2::ggplot(
  data = NULL,
  mapping = ggplot2::aes(
    x = time,
    y = intensity,
    group = mz,
    # color = as.character(inchikey_2D)
    color = molecular_formula
  )
) +
  ggplot2::geom_line(data = cc_pos) +
  ggplot2::geom_line(data = cc_neg) +
  ggplot2::geom_line(
    data = cc_cad_pos,
    color = "black",
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    data = cc_cad_neg,
    color = "black",
    linetype = "dashed"
  ) +
  ggplot2::scale_color_manual(
    values = c("#4E79A7", "#BAB0AC"),
    name = "Molecular Formula"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    complete = FALSE,
    # legend.position = "none",
    validate = TRUE
  ) +
  xlab("Time [min]") +
  ylab("Normalized Intensity")

plot_ik <- ggplot2::ggplot(
  data = NULL,
  mapping = ggplot2::aes(
    x = time,
    y = intensity,
    group = mz,
    # color = as.character(inchikey_2D)
    color = inchikey_2D
  )
) +
  ggplot2::geom_line(data = cc_pos |>
    dplyr::mutate(inchikey_2D = gsub("MIJYXULNPSFWEK", "oleanolic acid", inchikey_2D)) |>
    dplyr::mutate(inchikey_2D = gsub("WCGUUGGRBIKTOS", "ursolic acid", inchikey_2D))) +
  ggplot2::geom_line(data = cc_neg |>
    dplyr::mutate(inchikey_2D = gsub("MIJYXULNPSFWEK", "oleanolic acid", inchikey_2D)) |>
    dplyr::mutate(inchikey_2D = gsub("WCGUUGGRBIKTOS", "ursolic acid", inchikey_2D))) +
  ggplot2::geom_line(
    data = cc_cad_pos,
    color = "black",
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    data = cc_cad_neg,
    color = "black",
    linetype = "dashed"
  ) +
  ggplot2::scale_color_manual(
    values = c("#4E79A7", "#BAB0AC", "#A0CBE8"),
    name = "Structure (no stereochemistry)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    complete = FALSE,
    # legend.position = "none",
    validate = TRUE
  ) +
  xlab(label = "Time [min]") +
  ylab(label = "Normalized Intensity")

#' export
ggplot2::ggsave(
  filename = "~/git/cascade/data/paper/cascade-4.pdf",
  plot = ggpubr::ggarrange(
    plot_comparison,
    plot_taxo,
    plot_mf,
    plot_ik,
    labels = "AUTO",
    align = "hv"
  ),
  width = 9,
  height = 9,
  limitsize = FALSE
)
# ggplot2::ggsave(filename = "~/git/cascade/data/paper/cascade-6.pdf",
#                 plot = fig_taxo)
# ggplot2::ggsave(
#   filename = "~/git/cascade/data/paper/cascade-5.pdf",
#   plot = ggpubr::ggarrange(# plots_cad_pos[[1]],
#     plots_cad_pos[[3]],
#     # labels = "AUTO",
#     # nrow = 2,
#     common.legend = FALSE,
#     legend = "bottom"),
#   width = 8,
#   height = 4.5,
#   limitsize = FALSE
# )
# ggplot2::ggsave(
#   filename = "~/git/cascade/data/paper/cascade-7-a.pdf",
#   plot = fig_minmaj,
#   width = 12,
#   height = 6,
#   limitsize = FALSE
# )
# plotly::save_image(
#   p = treemaps$special,
#   file = "data/paper/cascade-7-b.pdf",
#   width = 1600,
#   height = 900
# )
# readr::write_tsv(x = plots_2_cad[["table_taxo_maj_cor_conf_signal_1"]],
#                  file = "data/paper/final_table.tsv")
# readr::write_csv(
#   x = candidates_confident |>
#     dplyr::filter(mode == "pos"),
#   file = "data/paper/confident_candidates.csv",
#   na = ""
# )

end <- Sys.time()
message("Script finished in", format(end - start))
