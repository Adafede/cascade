#' Save treemaps progress
#'
#' @param xs XS
#' @param type Type
#'
#' @return Saved treemaps
#'
#' @examples NULL
save_treemaps_progress <- function(xs, type = "treemap") {
  p <- progressr::progressor(along = xs)
  setNames(
    object = xs,
    nm = xs
  ) |>
    furrr::future_map(
      .f = function(x) {
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
