#' Title
#'
#' @param xs
#'
#' @return
#' @export
#'
#' @examples
save_histograms_progress <- function(xs) {
  p <- progressr::progressor(along = xs)
  future.apply::future_lapply(
    X = setNames(object = xs,
                 nm = xs),
    FUN = function(x) {
      ggplot2::ggsave(
        plot = histograms[[x]],
        filename = file.path(
          export_dir_histograms,
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
        height = 9 * max((size[[x]] / 100), 1) * max((dimensions[[x]] / 30), 1),
        limitsize = FALSE
      )
    }
  )
}