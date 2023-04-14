library(future)
library(progressr)
strat <- ifelse(test = .Platform$OS.type == "unix",
  yes = "multicore",
  no = "multisession"
)
plan(strategy = strat, workers = nbrOfWorkers())
# handlers(global = TRUE)
handlers(
  handler_txtprogressbar(enable = TRUE),
  handler_progress(format = ":spin [:bar] ETA: :eta :percent")
)
