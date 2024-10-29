#' Make confident
#'
#' @param df Dataframe
#' @param score Score
#'
#' @return A dataframe containing annotations with scores above the confidence threshold set
#'
#' @examples NULL
make_confident <- function(df, score) {
  df_ready <- df |>
    dplyr::mutate(
      best_candidate_1 = ifelse(
        test = as.numeric(score_final) >= score,
        yes = best_candidate_1,
        no = "notConfident"
      ),
      best_candidate_2 = ifelse(
        test = as.numeric(score_final) >= score,
        yes = best_candidate_2,
        no = "notConfident notConfident"
      ),
      best_candidate_3 = ifelse(
        test = as.numeric(score_final) >= score,
        yes = best_candidate_3,
        no = "notConfident notConfident notConfident"
      )
    )

  return(df_ready)
}
