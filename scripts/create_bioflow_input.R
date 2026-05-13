ensure_cran_packages <- function(packages, repos = "https://cloud.r-project.org") {
  missing <- packages[!vapply(
    packages,
    function(pkg) requireNamespace(pkg, quietly = TRUE),
    FUN.VALUE = logical(1)
  )]

  if (length(missing) > 0) {
    install.packages(missing, repos = repos)
  }
}

read_table_auto <- function(path) {
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

  phenotype <- read_table_auto(phenotype_file)
  pedigree <- read_table_auto(pedigree_file)
  vcf <- vcfR::read.vcfR(genotype_vcf_file, verbose = FALSE)
  genotype <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)

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
  run_from_cli()
}
