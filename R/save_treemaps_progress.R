#' Save treemaps progress
#'
#' @param xs XS
#' @param type Type
#'
#' @return Saved treemaps
#'
#' @examples NULL
save_treemaps_progress <- function(xs, type = "treemap") {
  setNames(
    object = xs,
    nm = xs
  ) |>
    purrr::map(
      .progress = TRUE,
      .f = function(x) {
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
