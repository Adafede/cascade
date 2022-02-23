#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
tables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = xs,
    FUN = function(x) {
      p()
      if (nrow(x != 0)) {
        x %>%
          dplyr::left_join(structures_classified) %>%
          dplyr::mutate(structureImage = RCurl::curlEscape(structureSmiles)) %>%
          dplyr::relocate(structureImage, .after = structure) %>%
          dplyr::relocate(structureLabel, .before = structure) %>%
          dplyr::select(-references_ids, -structure_id, -structureSmiles) %>%
          splitstackshape::cSplit(
            c("taxa", "taxaLabels", "references", "referencesLabels"),
            sep = "|",
            direction = "long"
          ) %>%
          dplyr::group_by(structure) %>%
          tidyr::fill(c("taxa", "taxaLabels", "references", "referencesLabels"),
            .direction = "downup"
          ) %>%
          dplyr::group_by(chemical_class) %>%
          dplyr::add_count(sort = TRUE) %>%
          dplyr::select(-n) %>%
          dplyr::group_by(chemical_superclass) %>%
          dplyr::add_count(sort = TRUE) %>%
          dplyr::select(-n) %>%
          dplyr::group_by(chemical_pathway) %>%
          dplyr::add_count(sort = TRUE) %>%
          dplyr::select(-n) %>%
          dplyr::distinct()
      } else {
        data.frame() |>
          dplyr::mutate(
            structureLabel = NA,
            structure = NA,
            structureImage = NA,
            taxaLabels = NA,
            taxa = NA,
            referenceLabels = NA,
            references = NA,
            chemical_pathway = NA,
            chemical_superclass = NA,
            chemical_class = NA
          )
      }
    }
  )
}
