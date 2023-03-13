library(future)
library(progressr)

plan(multisession)
# handlers(global = TRUE)
handlers(
  handler_txtprogressbar(enable = TRUE),
  handler_progress(format = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta")
)
