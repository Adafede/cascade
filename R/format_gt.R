require(package = dplyr, quietly = TRUE)
require(package = gt, quietly = TRUE)
require(package = htmltools, quietly = TRUE)
require(package = purrr, quietly = TRUE)

#' Title
#'
#' @param table
#'
#' @return
#' @export
#'
#' @examples
temp_gt_function <- function(table, title, subtitle) {
  prettyTable <- table %>%
    dplyr::mutate(
      structure = purrr::map(structure, ~ htmltools::a(href = .x, .x)),
      structure = purrr::map(structure, ~ gt::html(as.character(.x)))
    ) %>%
    dplyr::mutate(
      taxa = purrr::map(taxa, ~ htmltools::a(href = .x, .x)),
      taxa = purrr::map(taxa, ~ gt::html(as.character(.x)))
    ) %>%
    dplyr::mutate(
      references = purrr::map(references, ~ htmltools::a(href = .x, .x)),
      references = purrr::map(references, ~ gt::html(as.character(.x)))
    ) %>%
    gt::gt() %>%
    # cols_width(
    #   everything() ~ px(200)
    # ) %>%
    gt::tab_header(
      title = gt::md(title),
      subtitle = gt::md(subtitle)
    ) %>%
    gt::text_transform(
      locations = gt::cells_body(columns = structureImage),
      fn = molinfo
    )
  return(prettyTable)
}
