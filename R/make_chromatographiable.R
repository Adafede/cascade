#' Make chromatographiable
#'
#' @param df Dataframe
#' @param mass_min Mass min
#' @param mass_max Mass max
#' @param logp_min Log P min
#' @param logp_max Log P max
#'
#' @return A dataframe containing chromatographiable compounds
#'
#' @examples NULL
make_chromatographiable <-
  function(df,
           mass_min = 50,
           mass_max = 1500,
           logp_min = -1,
           logp_max = 6) {
    message("Keeping chromatographiable structures only")
    df_ready <- df |>
      dplyr::filter(structure_exact_mass >= mass_min) |>
      dplyr::filter(structure_exact_mass <= mass_max) |>
      dplyr::filter(structure_xlogp >= logp_min) |>
      dplyr::filter(structure_xlogp <= logp_max)

    return(df_ready)
  }
