bundle_bioflow_scripts <- function(
  input_files = c(
    "packages_verification.R",
    "read_geno_functions.R",
    "pheno_column_mapping_hardcoded.R",
    "pheno_column_mapping_utils.R",
    "getBioflowRdata.R"
  ),
  scripts_dir = "scripts",
  output_file = "bundled_getBioflowRdata.R",
  mapping_json_file = "resources/mapping.json",
  mapping_output_file = "pheno_column_mapping_hardcoded.R"
) {
  scripts_dir <- normalizePath(scripts_dir, mustWork = TRUE)

  if (!is.null(mapping_json_file) && nzchar(mapping_json_file)) {
    utils_script_path <- file.path(scripts_dir, "pheno_column_mapping_utils.R")
    if (!file.exists(utils_script_path)) {
      stop(
        sprintf("Missing mapping utility script: %s", utils_script_path),
        call. = FALSE
      )
    }

    source(utils_script_path, local = FALSE)
    write_mapping_fn <- get("write_hardcoded_column_mapping", mode = "function")
    write_mapping_fn(
      mapping_json_file = mapping_json_file,
      output_r_file = file.path(scripts_dir, mapping_output_file)
    )
  }

  input_paths <- file.path(scripts_dir, input_files)
  missing_files <- input_paths[!file.exists(input_paths)]
  if (length(missing_files) > 0) {
    stop(
      sprintf(
        "Missing input script(s): %s",
        paste(basename(missing_files), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  output_path <- file.path(scripts_dir, output_file)

  header <- c(
    "# Auto-generated file. Do not edit directly.",
    "# Source scripts are maintained in modular files under scripts/.",
    sprintf("# Generated on: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ""
  )

  blocks <- lapply(seq_along(input_paths), function(i) {
    block_title <- sprintf("# ---- BEGIN: %s ----", basename(input_paths[i]))
    block_footer <- sprintf("# ---- END: %s ----", basename(input_paths[i]))
    content <- readLines(input_paths[i], warn = FALSE, encoding = "UTF-8")
    c(block_title, content, block_footer, "")
  })

  bundled_content <- c(header, unlist(blocks, use.names = FALSE))
  writeLines(bundled_content, output_path, useBytes = TRUE)

  message(sprintf("Bundled script generated: %s", output_path))
  invisible(output_path)
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  scripts_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "scripts"
  output_file <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "bundled_getBioflowRdata.R"
  mapping_json_file <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL

  bundle_bioflow_scripts(
    scripts_dir = scripts_dir,
    output_file = output_file,
    mapping_json_file = mapping_json_file
  )
}
