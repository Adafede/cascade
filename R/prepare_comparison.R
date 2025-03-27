#' Prepare comparison
#'
#' @param features_informed Features informed
#' @param features_not_informed Features not informed
#' @param candidates_confident Candidates confident
#' @param min_similarity_prefilter Min similarity pre filter
#' @param min_similarity_filter Min similarity filter
#' @param mode Mode
#' @param show_example Show example? Default to FALSE
#'
#' @return A list of peaks
#'
#' @examples NULL
prepare_comparison <- function(
  features_informed = NULL,
  features_not_informed = NULL,
  candidates_confident,
  min_similarity_prefilter = 0.6,
  min_similarity_filter = 0.8,
  mode = "pos",
  show_example = FALSE
) {
  peaks_compared <- features_informed |>
    load_features_informed(show_example = show_example) |>
    tidytable::mutate(mode = mode)
  peaks_outside <- features_not_informed |>
    load_features_not_informed(show_example = show_example) |>
    tidytable::mutate(mode = mode) |>
    tidytable::mutate(tidytable::across(
      tidytable::all_of(
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
      ),
      as.numeric
    ))

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
    tidytable::filter(comparison_score >= min_similarity_prefilter) ## TODO check negative values

  peaks_min_precor <- peaks_maj |>
    tidytable::anti_join(peaks_maj_precor) |>
    temp_fix_2() |>
    tidytable::bind_rows(peaks_min)

  message("keeping multiple features only if none was reported in the species")
  peaks_maj_precor_taxo <- peaks_maj_precor |>
    tidytable::rowwise() |>
    tidytable::mutate(
      taxo = tidytable::if_else(
        condition = grepl(pattern = best_candidate_organism, x = species),
        true = 1,
        false = 0
      )
    ) |>
    tidytable::mutate(
      taxo = tidytable::if_else(
        condition = is.na(taxo),
        true = 0,
        false = taxo
      )
    ) |>
    tidytable::mutate(
      taxo_2 = tidytable::if_else(
        condition = best_candidate_organism %in% species,
        true = 1,
        false = 0
      )
    ) |>
    tidytable::group_by(sample, peak_id) |> ## TODO switch to ID if needed
    tidytable::mutate(sum = sum(taxo)) |>
    tidytable::mutate(sum_2 = sum(taxo_2)) |>
    tidytable::mutate(
      keep = tidytable::if_else(
        condition = taxo_2 == 1,
        true = "Y",
        false = tidytable::if_else(
          condition = taxo == 1,
          true = tidytable::if_else(
            condition = sum != sum_2 & sum_2 == 0,
            true = "Y",
            false = "N"
          ),
          false = tidytable::if_else(
            condition = sum == 0,
            true = "Y",
            false = "N"
          )
        )
      )
    ) |>
    tidytable::filter(keep == "Y") |> ## TODO decide if genus
    tidytable::ungroup()

  peaks_min_precor_taxo <- peaks_maj |>
    tidytable::anti_join(peaks_maj_precor_taxo) |>
    temp_fix_2() |>
    tidytable::bind_rows(peaks_min)

  message(
    "keeping peaks similarities with score above ",
    min_similarity_filter
  )
  peaks_maj_precor_taxo_cor <- peaks_maj_precor_taxo |>
    tidytable::filter(comparison_score >= min_similarity_filter)

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
