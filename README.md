# EBS-hackathon-26
Here we will store the code produced in the EBS-hackathon-26

## Bioflow input conversion script

Use `scripts/bundled_getBioflowRdata.R` to convert phenotype, pedigree, and genotype (VCF) files into a Bioflow-compatible R object.

### What it does
- Checks for required CRAN packages and installs missing ones (currently `vcfR`, `adegenet`, `cli`, `rlang`)
- Reads phenotype and pedigree tabular files (CSV)
- Reads genotype data from VCF it can be compressed (.gz) or plain text
- Creates a `bioflow_input` R object and writes it as an `.robj` file

### Usage


For use you can use the script located in `scripts/bundled_getBioflowRdata.R`

```R
library("https://raw.githubusercontent.com/Breeding-Analytics/EBS-hackathon-26/refs/heads/main/scripts/bundled_getBioflowRdata.R")

# Paths for phenotype, pedigree and genotype data
phenotype_file <- "test/pheno.csv"
pedigree_file <- "test/PedF1.csv"
genotype_vcf_file <- "test/inputF1.vcf"

# Paths of output .RData
output_file <- "bioflow_input_test" # Filename of the output file
output_path <- "test"

# List target traits to analyze present on phenotype dataset
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
```

## Column mapping and bundling

Create a JSON file with old-to-new phenotype column names:

```json
{
  "Year": "year",
  "Field_Location": "site",
  "Replicate": "rep"
}
```

Generate/update the hardcoded named list used by the scripts:

```bash
Rscript -e "source('scripts/pheno_column_mapping_utils.R'); write_hardcoded_column_mapping('scripts/pheno_mapping.json')"
```

Bundle all modular scripts into one file:

```bash
Rscript scripts/bundle_scripts.R
```

Regenerate mapping from JSON and bundle in one command:

```bash
Rscript scripts/bundle_scripts.R scripts bundled_getBioflowRdata.R scripts/pheno_mapping.json
```
