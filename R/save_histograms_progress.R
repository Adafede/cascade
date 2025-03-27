#' Save histograms progress
#'
#' @param xs XS
#'
#' @return Saved histograms
#'
#' @examples NULL
save_histograms_progress <- function(xs) {
  setNames(
    object = xs,
    nm = xs
  ) |>
    purrr::map(
      .progress = TRUE,
      .f = function(x) {
        ggplot2::ggsave(
          plot = histograms[[x]],
          filename = file.path(
            paths$data$histograms$path,
            paste0(
              "histogram_",
              gsub(
                pattern = " ",
                replacement = "_",
                x = x
              ),
              ".pdf"
            )
          ),
          width = 16 * max((dimensions[[x]] / 30), 1),
          height = 9 *
            max((size[[x]] / 100), 1) *
            max((dimensions[[x]] / 30), 1),
          limitsize = FALSE
        )
      }
    )
}
