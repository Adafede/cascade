#' Check and load packages 1
#'
#' @param cran CRAN
#' @param bioconductor Bioconductor
#' @param github Github
#'
#' @return Log of installed packages
#'
#' @examples NULL
check_and_load_packages_1 <- function(cran = packages_cran,
                                      bioconductor = packages_bioconductor,
                                      github = packages_github) {
  installed_packages <- rownames(installed.packages())

  if (!is.null(bioconductor)) {
    cran <- cran |>
      append("BiocManager")
  }
  if (!is.null(github)) {
    cran <- cran |>
      append("remotes")
  }

  installed_packages_cran <- cran %in% installed_packages

  return(lapply(X = cran[!installed_packages_cran], FUN = install.packages) |>
    invisible())
}

#' Check and load packages 2
#'
#' @param cran CRAN
#' @param bioconductor Bioconductor
#' @param github Github
#'
#' @return Log of installed packages
#'
#' @examples NULL
check_and_load_packages_2 <- function(cran = packages_cran,
                                      bioconductor = packages_bioconductor,
                                      github = packages_github) {
  installed_packages <- rownames(installed.packages())
  installed_packages_cran <- cran %in% installed_packages
  installed_packages_bioconductor <-
    bioconductor %in% installed_packages
  installed_packages_github <- github %in% installed_packages

  if (!is.null(bioconductor)) {
    cran <- cran |>
      append("BiocManager")
  }
  if (!is.null(github)) {
    cran <- cran |>
      append("remotes")
  }

  if (any(installed_packages_cran == FALSE)) {
    install.packages(cran[!installed_packages_cran])
  }
  if (any(installed_packages_bioconductor == FALSE)) {
    BiocManager::install(bioconductor[!installed_packages_bioconductor])
  }
  if (any(installed_packages_github == FALSE)) {
    lapply(X = github[!installed_packages_github], FUN = remotes::install_github)
  }

  return(lapply(
    c(
      cran,
      bioconductor,
      gsub(
        pattern = ".*/",
        replacement = "",
        x = github
      )
    ),
    require,
    character.only = TRUE
  ) |>
    invisible())
}
