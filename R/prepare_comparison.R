prepare_comparison <- function(detector = "cad") {
  log_debug(x = "loading compared peaks")
  peaks_compared <-
    readr::read_delim(file = file.path(EXPORT_DIR, switch(detector,
      "bpi" = EXPORT_FILE_BPI,
      "cad" = EXPORT_FILE_CAD,
      "pda" = EXPORT_FILE_PDA
    )))

  peaks_outside <-
    readr::read_delim(file = file.path(EXPORT_DIR, switch(detector,
      "bpi" = EXPORT_FILE_BPI_2,
      "cad" = EXPORT_FILE_CAD_2,
      "pda" = EXPORT_FILE_PDA_2
    )))

  peaks_all <- peaks_compared |>
    dplyr::bind_rows(peaks_outside)

  log_debug(x = "joining compared peaks and candidates")
  peaks_maj <- peaks_compared |>
    dplyr::left_join(candidates_confident)

  peaks_min <- peaks_outside |>
    dplyr::left_join(candidates_confident)

  log_debug(x = "temporary fix") #' TODO
  peaks_maj <- peaks_maj |>
    dplyr::mutate(
      id = sample,
      integral = peak_area,
      intensity = feature_area
    )
  peaks_min <- peaks_min |>
    dplyr::mutate(
      id = sample,
      integral = peak_area,
      intensity = feature_area,
      peak_rt_apex = rt,
      peak_area = 1
    )

  log_debug(x = "keeping peaks similarities above desired (pre-)threshold only")
  peaks_maj_precor <- peaks_maj |>
    dplyr::filter(comparison_score >= PEAK_SIMILARITY_PREFILTER) #' TODO check negative values

  log_debug(x = "keeping multiple features only if none was reported in the species")
  peaks_maj_precor_taxo <- peaks_maj_precor |>
    dplyr::rowwise() |>
    dplyr::mutate(taxo = ifelse(
      test = grepl(
        pattern = best_candidate_organism,
        x = species
      ),
      yes = 1,
      no = 0
    )) |>
    dplyr::mutate(taxo = ifelse(
      test = is.na(taxo),
      yes = 0,
      no = taxo
    )) |>
    dplyr::mutate(taxo_2 = ifelse(
      test = best_candidate_organism %in% species,
      yes = 1,
      no = 0
    )) |>
    dplyr::group_by(sample, peak_id) |> #' TODO switch to ID if needed
    dplyr::mutate(sum = sum(taxo)) |>
    dplyr::mutate(sum_2 = sum(taxo_2)) |>
    dplyr::mutate(keep = ifelse(
      test = taxo_2 == 1,
      yes = "Y",
      no = ifelse(
        test = taxo == 1,
        yes = ifelse(
          test = sum != sum_2 & sum_2 == 0,
          yes = "Y",
          no = "N"
        ),
        no = ifelse(
          test = sum == 0,
          yes = "Y",
          no = "N"
        )
      )
    )) |>
    dplyr::filter(keep == "Y") |> #' TODO decide if genus
    dplyr::ungroup()

  log_debug(x = "keeping peaks similarities with score above", PEAK_SIMILARITY)
  peaks_maj_precor_taxo_cor <- peaks_maj_precor_taxo |>
    dplyr::filter(comparison_score >= PEAK_SIMILARITY)

  returned_list <- list(
    peaks_all,
    peaks_maj,
    peaks_min,
    peaks_maj_precor,
    peaks_maj_precor_taxo,
    peaks_maj_precor_taxo_cor
  )
  names(returned_list) <- c(
    "peaks_all",
    "peaks_maj",
    "peaks_min",
    "peaks_maj_precor",
    "peaks_maj_precor_taxo",
    "peaks_maj_precor_taxo_cor"
  )

  return(returned_list)
}
