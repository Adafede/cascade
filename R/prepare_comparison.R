prepare_comparison <- function(detector = "cad") {
  log_debug(x = "loading compared peaks")
  path_1 <- switch(detector,
    "bpi" = IMPORT_FILE_BPI,
    "cad" = IMPORT_FILE_CAD,
    "pda" = IMPORT_FILE_PDA
  )
  path_2 <- switch(detector,
    "bpi" = IMPORT_FILE_BPI_2,
    "cad" = IMPORT_FILE_CAD_2,
    "pda" = IMPORT_FILE_PDA_2
  )
  peaks_compared <- path_1 |>
    lapply(
      FUN = function(x) {
        readr::read_delim(file = file.path(EXPORT_DIR, x)) |>
          # dplyr::mutate(peak_area = peak_area / max(peak_area))|>
          dplyr::mutate(mode = ifelse(
            test = grepl(
              pattern = "_pos_",
              x = x,
              ignore.case = TRUE
            ),
            yes = "pos",
            no = "neg"
          ))
      }
    ) |>
    dplyr::bind_rows()

  peaks_outside <- path_2 |>
    lapply(
      FUN = function(x) {
        readr::read_delim(file = file.path(EXPORT_DIR, x)) |>
          dplyr::mutate(mode = ifelse(
            test = grepl(
              pattern = "_pos_",
              x = x,
              ignore.case = TRUE
            ),
            yes = "pos",
            no = "neg"
          ))
      }
    ) |>
    dplyr::bind_rows()

  peaks_all <- peaks_compared |>
    dplyr::bind_rows(peaks_outside)

  log_debug(x = "joining compared peaks and candidates")
  log_debug(x = "temporary fix") ## TODO
  temp_fix <- function(df) {
    df_temp <- df |>
      dplyr::left_join(candidates_confident) |>
      dplyr::mutate(rt = as.numeric(rt))
    return(df_temp)
  }
  temp_fix_2 <- function(df) {
    df_temp <- df |>
      dplyr::mutate(
        id = sample,
        integral = peak_area,
        intensity = feature_area
      )
    return(df_temp)
  }
  temp_fix_3 <- function(df) {
    df_temp <- df |>
      dplyr::mutate(
        peak_rt_apex = rt,
        peak_area = 0.001
      )
    return(df_temp)
  }

  peaks_maj <- peaks_compared |>
    temp_fix() |>
    temp_fix_2()

  # peaks_maj_2 <- peaks_maj |>
  #   dplyr::rowwise() |>
  #   dplyr::mutate(acn = 5 + 95 / (acn_time / (pmin(
  #     acn_time, (peak_rt_apex - dvol - 0.5)
  #   )))) %>%
  #   dplyr::mutate(peak_area_corrected = predict_response(acn = acn, peak_area = peak_area))
  #
  # test <- peaks_maj_2 |>
  #   distinct(peak_id, peak_area, peak_rt_apex, peak_area_corrected)
  #
  # plotly::plot_ly(
  #   test,
  #   x =  ~ peak_rt_apex,
  #   y =  ~ peak_area,
  #   type = "bar"
  # ) |>
  #   plotly::add_bars(x =  ~ peak_rt_apex, y =  ~ peak_area_corrected)

  peaks_min <- peaks_outside |>
    temp_fix() |>
    temp_fix_2() |>
    temp_fix_3()

  log_debug(x = "keeping peaks similarities above desired (pre-)threshold only")
  peaks_maj_precor <- peaks_maj |>
    dplyr::filter(comparison_score >= PEAK_SIMILARITY_PREFILTER) ## TODO check negative values

  peaks_min_precor <- peaks_maj |>
    dplyr::anti_join(peaks_maj_precor) |>
    temp_fix_2() |>
    dplyr::bind_rows(peaks_min)

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
    dplyr::group_by(sample, peak_id) |> ## TODO switch to ID if needed
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
    dplyr::filter(keep == "Y") |> ## TODO decide if genus
    dplyr::ungroup()

  peaks_min_precor_taxo <- peaks_maj |>
    dplyr::anti_join(peaks_maj_precor_taxo) |>
    temp_fix_2() |>
    dplyr::bind_rows(peaks_min)

  log_debug(x = "keeping peaks similarities with score above", PEAK_SIMILARITY)
  peaks_maj_precor_taxo_cor <- peaks_maj_precor_taxo |>
    dplyr::filter(comparison_score >= PEAK_SIMILARITY)

  peaks_min_precor_taxo_cor <- peaks_maj |>
    dplyr::anti_join(peaks_maj_precor_taxo_cor) |>
    temp_fix_2() |>
    dplyr::bind_rows(peaks_min)

  returned_list <- list(
    peaks_all,
    peaks_maj,
    peaks_min,
    peaks_maj_precor,
    peaks_min_precor,
    peaks_maj_precor_taxo,
    peaks_min_precor_taxo,
    peaks_maj_precor_taxo_cor,
    peaks_min_precor_taxo_cor
  )
  names(returned_list) <- c(
    "peaks_all",
    "peaks_maj",
    "peaks_min",
    "peaks_maj_precor",
    "peaks_min_precor",
    "peaks_maj_precor_taxo",
    "peaks_min_precor_taxo",
    "peaks_maj_precor_taxo_cor",
    "peaks_min_precor_taxo_cor"
  )

  return(returned_list)
}
