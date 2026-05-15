load_json_column_mapping <- function(mapping_json_file) {
  if (!file.exists(mapping_json_file)) {
    stop(sprintf("Mapping JSON file not found: %s", mapping_json_file), call. = FALSE)
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read mapping JSON files.", call. = FALSE)
  }

  mapping_list <- jsonlite::fromJSON(mapping_json_file, simplifyVector = TRUE)

  if (!is.list(mapping_list) && !is.vector(mapping_list)) {
    stop("Invalid JSON mapping: expected an object with key/value pairs.", call. = FALSE)
  }

  mapping_vector <- unlist(mapping_list, use.names = TRUE)

  if (is.null(names(mapping_vector)) || any(names(mapping_vector) == "")) {
    stop("Invalid JSON mapping: all keys must be non-empty column names.", call. = FALSE)
  }

  if (!is.character(mapping_vector)) {
    mapping_vector <- as.character(mapping_vector)
  }

  if (anyNA(mapping_vector) || any(mapping_vector == "")) {
    stop("Invalid JSON mapping: all values must be non-empty target column names.", call. = FALSE)
  }

  mapping_vector
}

write_hardcoded_column_mapping <- function(
  mapping_json_file,
  output_r_file = file.path("scripts", "pheno_column_mapping_hardcoded.R"),
  object_name = "PHENO_COLUMN_MAPPING"
) {
  mapping_vector <- load_json_column_mapping(mapping_json_file)

  escaped_old <- gsub('"', '\\\\"', names(mapping_vector), fixed = TRUE)
  escaped_new <- gsub('"', '\\\\"', mapping_vector, fixed = TRUE)

  mapping_lines <- sprintf('  "%s" = "%s"', escaped_old, escaped_new)

  script_lines <- c(
    "# Auto-generated file. Do not edit directly.",
    sprintf("# Generated from JSON mapping: %s", normalizePath(mapping_json_file, winslash = "/", mustWork = TRUE)),
    "",
    sprintf("%s <- c(", object_name),
    paste(mapping_lines, collapse = ",\n"),
    ")",
    ""
  )

  writeLines(script_lines, output_r_file, useBytes = TRUE)
  message(sprintf("Hardcoded mapping script generated: %s", output_r_file))
  invisible(output_r_file)
}

apply_column_mapping_to_pheno <- function(data_pheno, column_mapping, strict = FALSE) {
  if (!is.data.frame(data_pheno)) {
    stop("data_pheno must be a data.frame.", call. = FALSE)
  }

  if (length(column_mapping) == 0) {
    return(data_pheno)
  }

  mapping_vector <- unlist(column_mapping, use.names = TRUE)

  if (is.null(names(mapping_vector)) || any(names(mapping_vector) == "")) {
    stop("column_mapping must be a named vector or named list.", call. = FALSE)
  }

  source_columns <- names(mapping_vector)
  target_columns <- as.character(mapping_vector)

  missing_in_input <- setdiff(source_columns, names(data_pheno))
  if (strict && length(missing_in_input) > 0) {
    stop(
      sprintf(
        "Mapped source column(s) missing in data_pheno: %s",
        paste(missing_in_input, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  rename_index <- match(names(data_pheno), source_columns)
  to_rename <- !is.na(rename_index)

  if (any(to_rename)) {
    names(data_pheno)[to_rename] <- target_columns[rename_index[to_rename]]
  }

  data_pheno
}
