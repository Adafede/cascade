#' Save pretty tables progress
#'
#' @param xs XS
#'
#' @return Saved pretty tables
#'
#' @examples NULL
save_prettyTables_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  setNames(
    object = xs,
    nm = xs
  ) |>
    furrr::future_map(
      .f = function(x) {
        p()
        gt::gtsave(
          data = prettyTables[[x]],
          filename = file.path(paths$data$tables$path, paste0(
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
