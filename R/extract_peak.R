#' Title
#'
#' @param peak
#'
#' @return
#' @export
#'
#' @examples
extract_peak <- function(peak) {
  rtr <- df_peak[peak, ] |>
    dplyr::mutate(
      rtmin = (rt_min + CAD_SHIFT) * 60,
      rtmax = (rt_max + CAD_SHIFT) * 60
    ) |>
    dplyr::select(rtmin, rtmax) |>
    as.matrix()

  mzr <- df_peak[peak, ] |>
    dplyr::select(mzmin = mz_min, mzmax = mz_max) |>
    as.matrix()

  ms_chr <- dda_data_min |>
    MSnbase::filterRt(c(rtr[1], rtr[2])) |>
    MSnbase::filterMz(c(mzr[1], mzr[2])) |>
    MSnbase::chromatogram(
      rt = rtr,
      mz = mzr
    )

  ms_ready <-
    data.frame(
      intensity = ms_chr[1]@intensity,
      rtime = ms_chr[1]@rtime
    ) |>
    dplyr::filter(!is.na(intensity)) |>
    dplyr::mutate(intensity = (intensity - min(intensity)) / (max(intensity) -
      min(intensity))) |>
    dplyr::filter(intensity >= 0.1) |>
    dplyr::mutate(rtime = (rtime - min(rtime)) / (max(rtime) -
      min(rtime))) |> #' see https://github.com/sneumann/xcms/issues/593
    dplyr::arrange(rtime)

  ms_peak <-
    MSnbase::Chromatogram(
      intensity = ms_ready$intensity,
      rtime = ms_ready$rtime
    )

  return(ms_peak)
}
