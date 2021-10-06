#' Title
#'
#' @param table
#'
#' @return
#' @export
#'
#' @examples
temp_gt_function <- function(table) {
  prettyTable <- table %>%
    mutate(
      structure = map(structure, ~ htmltools::a(href = .x, .x)),
      structure = map(structure, ~ gt::html(as.character(.x)))
    ) %>%
    mutate(
      taxa = map(taxa, ~ htmltools::a(href = .x, .x)),
      taxa = map(taxa, ~ gt::html(as.character(.x)))
    ) %>%
    mutate(
      references = map(references, ~ htmltools::a(href = .x, .x)),
      references = map(references, ~ gt::html(as.character(.x)))
    ) %>%
    gt() %>%
    # cols_width(
    #   everything() ~ px(200)
    # ) %>%
    text_transform(
      locations = cells_body(columns = structureImage),
      fn = molinfo
    )
  return(prettyTable)
}
