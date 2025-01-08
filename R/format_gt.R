#' Temp GT function
#'
#' @param table Table
#' @param title Title
#' @param subtitle Subtitle
#'
#' @return A formatted GT table
#'
#' @examples NULL
format_gt <- function(table,
                      title = "",
                      subtitle = "") {
  pretty_table <- table |>
    tidytable::mutate(
      structure = tidytable::map(structure, ~ htmltools::a(href = .x, .x)),
      structure = tidytable::map(structure, ~ gt::html(as.character(.x)))
    ) |>
    tidytable::mutate(
      taxa = tidytable::map(taxa, ~ htmltools::a(href = .x, .x)),
      taxa = tidytable::map(taxa, ~ gt::html(as.character(.x)))
    ) |>
    tidytable::mutate(
      references = tidytable::map(references, ~ htmltools::a(href = .x, .x)),
      references = tidytable::map(references, ~ gt::html(as.character(.x)))
    ) |>
    data.frame() |>
    gt::gt() |>
    # cols_width(
    #   everything() ~ px(200)
    # ) |>
    gt::tab_header(
      title = gt::md(title),
      subtitle = gt::md(subtitle)
    ) |>
    gt::text_transform(
      locations = gt::cells_body(columns = structureImage),
      fn = molinfo
    )
  return(pretty_table)
}
