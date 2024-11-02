#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
.datatable.aware <- TRUE

.onLoad <- function() {
  ## Conflicts with chromatographR
  # options("pboptions" = list(
  #   type = if (interactive()) "timer" else "none",
  #   char = "-",
  #   txt.width = 50,
  #   gui.width = 300,
  #   style = 3,
  #   initial = 0,
  #   title = "R progress bar",
  #   label = "",
  #   nout = 100L
  # ))
  strat <- ifelse(test = .Platform$OS.type == "unix",
    yes = "multicore",
    no = "multisession"
  )
  future::plan(strategy = strat, workers = future::nbrOfWorkers())
  progressr::handlers(
    progressr::handler_txtprogressbar(enable = TRUE),
    progressr::handler_progress(format = ":spin [:bar] ETA: :eta :percent")
  )
  invisible(NULL)
}

.onAttach <- function() {
  packageStartupMessage("Welcome to CASCADE")
}
