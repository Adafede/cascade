p1q1 <- 0.00001
p1q2 <- -0.0006
p1q3 <- -0.0778
p2q1 <- 0.00002
p2q2 <- -0.00022
p2q3 <- 0.05499
p3q1 <- -0.00017
p3q2 <- 0.0209
p3q3 <- 1.4041
acn <- 100
acn_time <- 30.5
dvol <- 0.5 ## dirty at the moment

#' Title
#'
#' @param acn_eluent
#' @param q1
#' @param q2
#' @param q3
#'
#' @return
#' @export
#'
#' @examples
p_acn_i <- function(acn_eluent, q1, q2, q3) {
  result <- q1 * acn_eluent^2 + q2 * acn_eluent + q3
  return(result)
}

#' Title
#'
#' @param acn
#' @param peak_area
#'
#' @return
#' @export
#'
#' @examples
predict_response <- function(acn, peak_area) {
  p_acn_1 <- p_acn_i(
    acn_eluent = acn,
    q1 = p1q1,
    q2 = p1q2,
    q3 = p1q3
  )

  p_acn_2 <- p_acn_i(
    acn_eluent = acn,
    q1 = p2q1,
    q2 = p2q2,
    q3 = p2q3
  )

  p_acn_3 <- p_acn_i(
    acn_eluent = acn,
    q1 = p3q1,
    q2 = p3q2,
    q3 = p3q3
  )

  concentration <- if (peak_area > 0) {
    peak_area^(1 / (2 * p_acn_1 + p_acn_2)) *
      exp(1)^(-p_acn_3 / (2 * p_acn_1 + p_acn_2))
  } else {
    0
  }

  return(concentration)
}
