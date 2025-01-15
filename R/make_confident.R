#' Make confident
#'
#' @param df Dataframe
#' @param score Score
#'
#' @return A dataframe containing annotations with scores above the confidence threshold set
#'
#' @examples NULL
make_confident <- function(df, score) {
  df |>
    tidytable::mutate(
      best_candidate_1 = tidytable::if_else(
        condition = as.numeric(score_final) >= score,
        true = best_candidate_1,
        false = "notConfident"
      ),
      best_candidate_2 = tidytable::if_else(
        condition = as.numeric(score_final) >= score,
        true = best_candidate_2,
        false = "notConfident notConfident"
      ),
      best_candidate_3 = tidytable::if_else(
        condition = as.numeric(score_final) >= score,
        true = best_candidate_3,
        false = "notConfident notConfident notConfident"
      )
    )
}
