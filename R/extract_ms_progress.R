#' Extract MS progress
#'
#' @param xs XS
#' @param ms_data MS Data
#' @param rts RTs
#' @param mzs MZs
#' @param nrows N rows
#'
#' @return A list of extracted MS peaks
#'
#' @examples NULL
extract_ms_progress <- function(xs, ms_data, rts, mzs, nrows) {
  safe_chromatogram <- function(ms_data, rt, mz, max_attempts = 10) {
    for (attempt in seq_len(max_attempts)) {
      tryCatch(
        expr = {
          return(
            MSnbase::chromatogram(
              ms_data,
              rt = rt,
              mz = mz,
              BPPARAM = BiocParallel::SerialParam()
            )
          )
        },
        error = function(e) {
          if (attempt < max_attempts) {
            warning(
              sprintf(
                "Chromatogram extraction failed (Attempt %d). Retrying in %d seconds. Error: %s",
                attempt,
                attempt,
                conditionMessage(e)
              )
            )
          } else {
            warning(
              sprintf(
                "Chromatogram extraction failed after %d attempts. Skipping. Error: %s",
                max_attempts,
                conditionMessage(e)
              )
            )
            return(NULL)
          }
        }
      )
    }
  }
  xs |>
    purrr::map(
      .progress = TRUE,
      .f = function(x, ms_data, rts, mzs, nrows) {
        message("CAD Peak: ", x)
        safe_chromatogram(
          ms_data = ms_data,
          rt = rts[[x]],
          mz = Reduce(f = rbind, x = mzs[[x]])
        ) |>
          transform_ms()
      },
      ms_data = ms_data,
      rts = rts,
      mzs = mzs,
      nrows = nrows
    )
}
