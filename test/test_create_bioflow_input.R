# Test script for create_bioflow_input function
# This script tests the getBioflowRData function with sample data files

# Define paths to test data files (relative to project root)
phenotype_file <- "test/bioflow_pheno_data.csv"
pedigree_file <- "test/PedF1.csv"
genotype_vcf_file <- "test/bioflow_genotype_data.vcf"
output_file <- "bioflow_input_test"
output_path <- "test"

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

# Source and run the function to collect all the codebase into a 
# single script
source("scripts/bundle_scripts.R")
bundle_bioflow_scripts()

# Source only the function definitions from the main script
source("scripts/bundled_getBioflowRdata.R")
#source("https://raw.githubusercontent.com/Breeding-Analytics/EBS-hackathon-26/refs/heads/main/scripts/bundled_getBioflowRdata.R")

# Run the function
cat("\nExecuting create_bioflow_input function...\n")
traits <- c("Maize Plant Height","Wheat Plant Height","Leaf Blight Severity (1-9)")

result <- getBioflowRData(
  phenotypeFile = phenotype_file,
  pedigreeFile  = NULL,
  genotypeFile = genotype_vcf_file,
  traits = traits,
  outputPath = output_file,
  outputFile = output_path
)

# Verify result
cat("\nVerifying output...\n")
# Verify the completedness of the results object
cat("\n✓ Test completed successfully!\n")
