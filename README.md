# EBS-hackathon-26
Here we will store the code produced in the EBS-hackathon-26

## Bioflow input conversion script

Use `scripts/create_bioflow_input.R` to convert phenotype, pedigree, and genotype (VCF) files into a Bioflow-compatible R object.

### What it does
- Checks for required CRAN packages and installs missing ones (currently `vcfR`)
- Reads phenotype and pedigree tabular files (CSV or TSV)
- Reads genotype data from VCF
- Creates a `bioflow_input` R object and writes it as an `.rds` file

### Usage

```bash
Rscript scripts/create_bioflow_input.R \
  /path/to/phenotype.csv \
  /path/to/pedigree.csv \
  /path/to/genotypes.vcf \
  /path/to/output/bioflow_input.rds
```
