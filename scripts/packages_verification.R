ensure_cran_packages <- function(packages, repos = "https://cloud.r-project.org") {
  missing <- packages[!vapply(
    packages,
    function(pkg) requireNamespace(pkg, quietly = TRUE),
    FUN.VALUE = logical(1)
  )]

  if (length(missing) > 0) {
    tryCatch(
      install.packages(missing, repos = repos),
      error = function(e) {
        stop(
          sprintf("Failed to install required CRAN packages: %s", paste(missing, collapse = ", ")),
          call. = FALSE
        )
      }
    )
  }

  still_missing <- packages[!vapply(
    packages,
    function(pkg) requireNamespace(pkg, quietly = TRUE),
    FUN.VALUE = logical(1)
  )]

  if (length(still_missing) > 0) {
    stop(
      sprintf("Required CRAN packages are still unavailable: %s", paste(still_missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

ensure_github_packages <- function(){
  if (!require(cgiarPipeline)) remotes::install_github("Breeding-Analytics/cgiarPipeline")
  if (!require(cgiarBase)) remotes::install_github("Breeding-Analytics/cgiarBase")
}