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

read_table_auto <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  }

  if (isTRUE(file.info(path)$size == 0)) {
    stop(sprintf("Input file is empty: %s", path), call. = FALSE)
  }

  first_line <- readLines(path, n = 1, warn = FALSE)
  sep <- if (grepl("\\t", first_line)) "\t" else ","

  utils::read.table(
    file = path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "\""
  )
}

create_bioflow_input <- function(phenotype_file,
                                 pedigree_file,
                                 genotype_vcf_file,
                                 output_file = NULL) {
  ensure_cran_packages(c("vcfR"))

  with_file_context <- function(label, path, reader) {
    tryCatch(
      reader(path),
      error = function(e) {
        stop(sprintf("Failed to read %s file '%s': %s", label, path, conditionMessage(e)), call. = FALSE)
      }
    )
  }

  phenotype <- with_file_context("phenotype", phenotype_file, read_table_auto)
  pedigree <- with_file_context("pedigree", pedigree_file, read_table_auto)
  vcf <- with_file_context("genotype VCF", genotype_vcf_file, function(path) {
    vcfR::read.vcfR(path, verbose = FALSE)
  })
  genotype <- tryCatch(
    vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE),
    error = function(e) {
      stop(
        sprintf("Failed to extract genotype matrix from VCF '%s': %s", genotype_vcf_file, conditionMessage(e)),
        call. = FALSE
      )
    }
  )

  bioflow_input <- list(
    phenotype = phenotype,
    pedigree = pedigree,
    genotype = genotype,
    metadata = list(
      phenotype_file = normalizePath(phenotype_file, mustWork = FALSE),
      pedigree_file = normalizePath(pedigree_file, mustWork = FALSE),
      genotype_vcf_file = normalizePath(genotype_vcf_file, mustWork = FALSE),
      created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  )

  class(bioflow_input) <- c("bioflow_input", class(bioflow_input))

  if (!is.null(output_file)) {
    saveRDS(bioflow_input, file = output_file)
  }

  bioflow_input
}

run_from_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4) {
    stop(
      paste(
        "Usage: Rscript scripts/create_bioflow_input.R",
        "<phenotype_file> <pedigree_file> <genotype_vcf_file> <output_rds_file>",
        sep = " "
      ),
      call. = FALSE
    )
  }

  create_bioflow_input(
    phenotype_file = args[[1]],
    pedigree_file = args[[2]],
    genotype_vcf_file = args[[3]],
    output_file = args[[4]]
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  tryCatch(
    run_from_cli(),
    error = function(e) {
      message(conditionMessage(e))
      quit(status = 1)
    }
  )
}
