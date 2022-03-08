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
           mass_min = params$structures$mass$min,
           mass_max = params$structures$mass$max,
           logp_min = params$structures$logp$min,
           logp_max = params$structures$logp$max) {
    message("Keeping chromatographiable structures only")
    df_ready <- df |>
      dplyr::filter(structure_exact_mass >= mass_min) |>
      dplyr::filter(structure_exact_mass <= mass_max) |>
      dplyr::filter(structure_xlogp >= logp_min) |>
      dplyr::filter(structure_xlogp <= logp_max)

    return(df_ready)
  }
