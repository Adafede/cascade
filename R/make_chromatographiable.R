#' Title
#'
#' @param df
#' @param mass_min
#' @param mass_max
#' @param logp_min
#' @param logp_max
#'
#' @return
#' @export
#'
#' @examples
make_chromatographiable <-
  function(df,
           mass_min = 100,
           mass_max = 1500,
           logp_min = -1,
           logp_max = 6) {
    df_ready <- df |>
      dplyr::filter(structure_exact_mass >= mass_min) |>
      dplyr::filter(structure_exact_mass <= mass_max) |>
      dplyr::filter(structure_xlogp >= logp_min) |>
      dplyr::filter(structure_xlogp <= logp_max)

    return(df_ready)
  }
