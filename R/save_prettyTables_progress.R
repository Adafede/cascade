#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
save_prettyTables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = setNames(
      object = xs,
      nm = xs
    ),
    FUN = function(x) {
      gt::gtsave(
        data = prettyTables[[x]],
        filename = file.path(export_dir_tables, paste0(
          "prettyTable_",
          gsub(
            pattern = " ",
            replacement = "_",
            x = x
          ),
          ".html"
        ))
      )
    }
  )
}
