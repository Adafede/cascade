#' Save pretty subtables progress
#'
#' @param xs XS
#'
#' @return Saved pretty subtables
#'
#' @examples NULL
save_prettySubtables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  setNames(
    object = xs,
    nm = xs
  ) |>
    furrr::future_map(
      .f = function(x) {
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
