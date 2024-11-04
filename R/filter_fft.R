#' Filter FFT
#'
#' @param x X
#' @param components Components
#'
#' @return The fourier filtered x
#'
#' @examples NULL
filter_fft <- function(x, components) {
  # nucleR::filterFFT(data = x, pcKeepComp = components, useOptim = TRUE)
  x[is.na(x)] <- 0
  temp <- fft(x)
  keep <- round(length(temp) * components)
  temp[keep:(length(temp) - keep)] <- 0
  return(Re(fft(temp, inverse = TRUE)) / length(temp))
}
