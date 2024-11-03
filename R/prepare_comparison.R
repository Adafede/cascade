#' Prepare comparison
#'
#' @param detector Detector
#'
#' @return A list of peaks
#'
#' @examples NULL
prepare_comparison <- function(detector = "cad") {
  message("loading compared peaks")
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
        tidytable::fread(file = file.path(x)) |>
          # tidytable::mutate(peak_area = peak_area / max(peak_area))|>
          tidytable::mutate(mode = ifelse(
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
    tidytable::bind_rows()

  peaks_outside <- path_2 |>
    lapply(
      FUN = function(x) {
        tidytable::fread(file = file.path(x)) |>
          tidytable::mutate(mode = ifelse(
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
    tidytable::bind_rows() |>
    tidytable::mutate(tidytable::across(tidytable::all_of(
      c(
        "peak_id",
        "peak_rt_min",
        "peak_rt_apex",
        "peak_rt_max",
        "peak_area",
        "feature_id",
        "feature_rt",
        "feature_mz",
        "feature_area",
        "comparison_score"
      )
    ), as.numeric))

  peaks_all <- peaks_compared |>
    tidytable::bind_rows(peaks_outside)

  message("joining compared peaks and candidates")
  message("temporary fix") ## TODO
  temp_fix <- function(df) {
    df |>
      tidytable::mutate(feature_id = as.numeric(feature_id)) |>
      tidytable::left_join(candidates_confident) |>
      tidytable::mutate(rt = as.numeric(rt))
  }

  temp_fix_2 <- function(df) {
    df |>
      tidytable::mutate(
        id = sample,
        integral = peak_area,
        intensity = feature_area
      )
  }

  temp_fix_3 <- function(df) {
    df |>
      tidytable::mutate(peak_rt_apex = rt, peak_area = 0.001)
  }

  temp_fix_4 <- function(df) {
    df |>
      tidytable::mutate(
        peak_id = as.numeric(peak_id),
        peak_rt_min = as.numeric(peak_rt_min),
        peak_rt_max = as.numeric(peak_rt_max),
        feature_rt = as.numeric(feature_rt),
        feature_mz = as.numeric(feature_mz),
        feature_area = as.numeric(feature_area),
        comparison_score = as.numeric(comparison_score),
        integral = as.numeric(integral),
        intensity = as.numeric(intensity)
      )
  }

  peaks_maj <- peaks_compared |>
    temp_fix() |>
    temp_fix_2()

  # peaks_maj_2 <- peaks_maj |>
  #   tidytable::rowwise() |>
  #   tidytable::mutate(acn = 5 + 95 / (acn_time / (pmin(
  #     acn_time, (peak_rt_apex - dvol - 0.5)
  #   )))) |>
  #   tidytable::mutate(peak_area_corrected = predict_response(acn = acn, peak_area = peak_area))
  #
  # test <- peaks_maj_2 |>
  #   tidytable::distinct(peak_id, peak_area, peak_rt_apex, peak_area_corrected)
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
    temp_fix_3() |>
    temp_fix_4()

  message("keeping peaks similarities above desired (pre-)threshold only")
  peaks_maj_precor <- peaks_maj |>
    tidytable::filter(comparison_score >= PEAK_SIMILARITY_PREFILTER) ## TODO check negative values

  peaks_min_precor <- peaks_maj |>
    tidytable::anti_join(peaks_maj_precor) |>
    temp_fix_2() |>
    tidytable::bind_rows(peaks_min)

  message("keeping multiple features only if none was reported in the species")
  peaks_maj_precor_taxo <- peaks_maj_precor |>
    tidytable::rowwise() |>
    tidytable::mutate(taxo = ifelse(
      test = grepl(pattern = best_candidate_organism, x = species),
      yes = 1,
      no = 0
    )) |>
    tidytable::mutate(taxo = ifelse(
      test = is.na(taxo),
      yes = 0,
      no = taxo
    )) |>
    tidytable::mutate(taxo_2 = ifelse(
      test = best_candidate_organism %in% species,
      yes = 1,
      no = 0
    )) |>
    tidytable::group_by(sample, peak_id) |> ## TODO switch to ID if needed
    tidytable::mutate(sum = sum(taxo)) |>
    tidytable::mutate(sum_2 = sum(taxo_2)) |>
    tidytable::mutate(keep = ifelse(
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
    tidytable::filter(keep == "Y") |> ## TODO decide if genus
    tidytable::ungroup()

  peaks_min_precor_taxo <- peaks_maj |>
    tidytable::anti_join(peaks_maj_precor_taxo) |>
    temp_fix_2() |>
    tidytable::bind_rows(peaks_min)

  message("keeping peaks similarities with score above", PEAK_SIMILARITY)
  peaks_maj_precor_taxo_cor <- peaks_maj_precor_taxo |>
    tidytable::filter(comparison_score >= PEAK_SIMILARITY)

  peaks_min_precor_taxo_cor <- peaks_maj |>
    tidytable::anti_join(peaks_maj_precor_taxo_cor) |>
    temp_fix_2() |>
    tidytable::bind_rows(peaks_min)

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
