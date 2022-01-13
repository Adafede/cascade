#' Title
#'
#' @param df
#' @param score
#'
#' @return
#' @export
#'
#' @examples
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
