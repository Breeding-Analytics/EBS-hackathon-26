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

'''
# Run the script on test/test_create_bioflow_input.R
'''
