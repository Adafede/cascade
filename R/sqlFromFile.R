#' SQL from file
#'
#' @source https://stackoverflow.com/questions/18914283/how-to-execute-more-than-one-rsqlite-statement-at-once-or-how-to-dump-a-whole-fi
#'
#' @param file File
#'
#' @return SQL todo
#'
#' @examples NULL
sqlFromFile <- function(file) {
  sql <- suppressWarnings(readLines(file))
  sql <- unlist(stringr::str_split(paste(sql, collapse = " "), ";"))
  sql <- sql[grep("^ *$", sql, invert = TRUE)]
  sql
}
