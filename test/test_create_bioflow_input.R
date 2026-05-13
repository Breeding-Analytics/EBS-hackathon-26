# Test script for create_bioflow_input function
# This script tests the create_bioflow_input function with sample data files

# Define paths to test data files
phenotype_file <- "test/pheno.csv"
pedigree_file <- "test/PedF1.csv"
genotype_vcf_file <- "test/inputF1.vcf"
output_file <- "test/bioflow_input_test.rds"

# Verify test data files exist
cat("Checking for test data files...\n")
if (!file.exists(phenotype_file)) {
  stop(sprintf("Phenotype file not found: %s", phenotype_file), call. = FALSE)
}
if (!file.exists(pedigree_file)) {
  stop(sprintf("Pedigree file not found: %s", pedigree_file), call. = FALSE)
}
if (!file.exists(genotype_vcf_file)) {
  stop(sprintf("Genotype VCF file not found: %s", genotype_vcf_file), call. = FALSE)
}
cat("✓ All test data files found\n")

# Source only the function definitions from the main script
source_lines <- readLines("scripts/create_bioflow_input.R")

# Remove the CLI execution section (last few lines)
# Keep everything before the "if (identical(environment()..." line
script_lines <- source_lines[1:(which(grepl("if \\(identical\\(environment", source_lines)) - 1)]
script_content <- paste(script_lines, collapse = "\n")

# Execute the functions
eval(parse(text = script_content))

# Run the function
cat("\nExecuting create_bioflow_input function...\n")
result <- create_bioflow_input(
  phenotype_file = phenotype_file,
  pedigree_file = pedigree_file,
  genotype_vcf_file = genotype_vcf_file,
  output_file = output_file
)

# Verify result
cat("\nVerifying output...\n")
cat(sprintf("✓ Result class: %s\n", paste(class(result), collapse = ", ")))
cat(sprintf("✓ Phenotype rows: %d\n", nrow(result$phenotype)))
cat(sprintf("✓ Pedigree rows: %d\n", nrow(result$pedigree)))
cat(sprintf("✓ Genotype dimensions: %d x %d\n", nrow(result$genotype), ncol(result$genotype)))
cat(sprintf("✓ Output saved to: %s\n", output_file))

cat("\n✓ Test completed successfully!\n")
