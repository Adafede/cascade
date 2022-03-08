#' Title
#'
#' @param xs
#' @param type
#'
#' @return
#' @export
#'
#' @examples
save_treemaps_progress <- function(xs, type = "treemap") {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = setNames(
      object = xs,
      nm = xs
    ),
    FUN = function(x) {
      p()
      plotly::save_image(
        p = switch(type,
          "treemap" = treemaps,
          "sunburst" = sunbursts
        )[[x]],
        file = file.path(
          switch(type,
            "treemap" = paths$data$treemaps$path,
            "sunburst" = paths$data$sunbursts$path
          ),
          paste0(
            type,
            "_",
            gsub(
              pattern = " ",
              replacement = "_",
              x = x
            ),
            ".pdf"
          )
        ),
        width = 900,
        height = 900
      )
    }
  )
}
