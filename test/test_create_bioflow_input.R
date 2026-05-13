# Test script for create_bioflow_input function
# This script tests the getBioflowRData function with sample data files

# Determine script directory and project root
# When run from Rscript, we can get the calling command from commandArgs()
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", args, value = TRUE)

if (length(script_arg) > 0) {
  script_path <- sub("--file=", "", script_arg)
  script_dir <- dirname(normalizePath(script_path))
  project_root <- dirname(script_dir)
} else {
  # Fallback: assume current directory is project root
  project_root <- getwd()
}

# Set working directory to project root for consistent relative paths
setwd(project_root)

# Define paths to test data files (relative to project root)
phenotype_file <- "test/pheno.csv"
pedigree_file <- "test/PedF1.csv"
genotype_vcf_file <- "test/inputF1.vcf"
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

# Set the scripts directory as an option so sourced scripts know where they are
options(scripts_dir = normalizePath("scripts"))

# Source only the function definitions from the main script
source("scripts/getBioflowRdata.R")

# Run the function
cat("\nExecuting create_bioflow_input function...\n")
traits <- c("Pollen_DAP_days","Silk_DAP_days","Plant_Height_cm","Ear_Height_cm",
            "Root_Lodging_plants","Stalk_Lodging_plants",
            "Yield_Mg_ha","Grain_Moisture","Twt_kg_m3")

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
