#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
save_prettySubtables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = setNames(
      object = xs,
      nm = xs
    ),
    FUN = function(x) {
      p()
      gt::gtsave(
        data = prettySubtables[[x]],
        filename = file.path(
          paths$data$tables$path,
          paste0(
            "prettySubtable_",
            gsub(
              pattern = " ",
              replacement = "_",
              x = x
            ),
            ".html"
          )
        )
      )
    }
  )
}
