# Copyright (C) 2026 Enterprise Breeding System
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# -------------------------------------------------------------------------------------
# Name             : getBioflowRData
# Description      : Create RData that will be use by Bioflow ready for Single Trial
#                    Analysis (STA)
# R Version        : 4.2.3
# Pkg Dependency   : 
# Note             : Revise how environment is created 
# ---------------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@cgiar.org
# Date             : 2024.10.15
# Date Modified    : 2026.05.13
# Maintainer       : Bioflow core development team 
# Maintainer Email : 
# Script Version   : 1.8
# Command          : getBioflowRData(phenotypeFile = "C:/sta_input_phenotype_v3_rev.csv"),
#                             traits = c("YLD_TON","FLW50"),
#                             outputPath = "C:/results",
#                             outputFile = "Request1",
#                             requestId = "1010")
# ---------------------------------------------------------------------------------------------
#' @name getBioflowRData
#' @aliases getBioflowRData
#' @title Create R object for Bioflow with phenotypic and genotypic QA/QC
#' 
#' @description Create RData object with Results structure with pheno and geno QA/QC ready for downstream analysis
#' @param phenotypeFile a string indicating the path and file name in csv format containing the QCed phenotypic data for analysis from EBS
#' @param pedigreeFile NULL or string indicating the path and file name in csv format containing the pedigree data associated with the phenotype file
#' @param genotypeFile NULL or string indicating the path and file name of VCF file containing the genotypic data associated with the phenotype file
#' @param traits a character vector indicating the trait that will be use for analysis
#' @param outputPath string indicating the path where the RData will be save
#' @param outputFile string indicating the RData file name that will be use 
#' @param requestId string indicating the analysis request id generated in EBS 
# -------------------------------------------------------------------------------------

source_relative <- function(filename) {
  # Get the scripts directory from options if set, otherwise try to determine it
  scripts_dir <- getOption("scripts_dir")
  
  if (is.null(scripts_dir)) {
    # Try to determine scripts directory from call stack
    for (i in 1:sys.nframe()) {
      frame <- sys.frame(i)
      srcfile <- attr(frame, "srcfile")
      if (!is.null(srcfile)) {
        file_path <- normalizePath(srcfile$filename)
        # If this file is in scripts directory, use that
        if (grepl("/scripts/", file_path) || grepl("\\\\scripts\\\\", file_path)) {
          scripts_dir <- dirname(file_path)
          break
        }
        # Otherwise check if parent directory has a scripts subdirectory
        parent_dir <- dirname(file_path)
        if (dir.exists(file.path(parent_dir, "scripts"))) {
          scripts_dir <- file.path(parent_dir, "scripts")
          break
        }
      }
    }
  }
  
  if (is.null(scripts_dir)) {
    stop("Could not determine scripts directory. Please set options(scripts_dir = '...')", call. = FALSE)
  }
  
  full_path <- file.path(scripts_dir, filename)
  if (!file.exists(full_path)) {
    stop(sprintf("File not found: %s (looked in: %s)", filename, scripts_dir), call. = FALSE)
  }
  source(full_path, local = FALSE)
}

getBioflowRData <- function(phenotypeFile, pedigreeFile = NULL, genotypeFile = NULL, 
                     traits, outputPath, outputFile, requestId = NULL) {
  
  # Source other required scripts using location-independent paths
  source_relative("packages_verification.R")
  
  # Now we can use ensure_cran_packages since packages_verification.R has been sourced
  ensure_cran_packages(c("vcfR", "adegenet", "cli", "rlang"))
  
  # Source read_geno_functions.R
  source_relative("read_geno_functions.R")
  
  analysisId <- round(as.numeric(Sys.time()), 0)
  cat(paste0("analysisId: ", analysisId, "\n"))
  
  traits <- unlist(traits)  # ensure object is a vector
  cat(paste0("trait: ", traits, "\n"))
  cat("\n")
  
  # --- data -> pheno
  # data_pheno <- read.csv(phenotypeFile, encoding = 'utf-8')
  
  # # remove non-alphanumeric character in occurrenceName and revise how enviroment is created
  # data_pheno$environment <- paste0("Env", paste(data_pheno$year, data_pheno$season, gsub("[^a-zA-Z0-9]", "", data_pheno$occurrenceName), sep = "_"))
  
  # # check if variable are present
  # if (any(is.na(data_pheno$design)) || any(data_pheno$design == "")) {
  #   stop("Design information are missing for some or all rows.")
  # }
  
  # # check design parameter base, rep in EBS pertains to the number of times that entry appear in the occurrence 
  # data_pheno[data_pheno$design == "Partially Replicated", "rep"] <- NA
  # data_pheno[data_pheno$design == "Augmented", "rep"] <- NA
  # data_pheno[data_pheno$design == "Augmented RCBD", "rep"] <- NA # assumes that blockNumber is not NA
  
  # # --- data -> pedigree
  # data_pedigree <- data.frame(
  #   germplasmName = as.vector(unique(data_pheno$germplasmName)),
  #   mother = NA,
  #   father = NA,
  #   yearOfOrigin = NA
  # )
  
  # --- data -> geno
  data_geno <- read_vcf(genotypeFile) # todo generalize the ploidity level
  # --- metadata -> pheno
  # NOTE: rep() is bugged with python
  # metadata_pheno <- data.frame(
  #   parameter = as.vector(c("stage", "year", "season", "location", "trial", "study", "rep", "iBlock", "row", "col", "designation", "gid", "entryType", rep("trait", nTrait), "environment", "source")),
  #   value = as.vector(c("breedingStage", "year", "season", "site", "experiment", "occurrenceName", "rep", "blockNumber", "paY", "paX", "germplasmName", "germplasmDbId", "entryType", traits, "environment", "ebs-ba"))
  # )
  
  # metadata_pheno1 <- data.frame(
  #   parameter = c("stage", "year", "season", "location", "trial", "study", "rep", "iBlock", "row", "col", "designation", "gid", "entryType"),
  #   value = c("breedingStage", "year", "season", "site", "experiment", "occurrenceName", "rep", "blockNumber", "paY", "paX", "germplasmName", "germplasmDbId", "entryType")
  # )
  # metadata_pheno2 <- data.frame(
  #   parameter = "trait",
  #   value = traits
  # )
  # if (is.null(requestId)) {
  #   metadata_pheno3 <- data.frame(
  #     parameter = c("environment", "source", "sourceId"),
  #     value = c("environment", "ebs-ba", NA)
  #   )
  # } else {
  #   metadata_pheno3 <- data.frame(
  #     parameter = c("environment", "source", "sourceId"),
  #     value = c("environment", "ebs-ba", requestId)
  #   )  
  # }
  
  # metadata_pheno <- rbind(rbind(metadata_pheno1, metadata_pheno2), metadata_pheno3)
  
  # metadata_pheno_parameter_size <- length(metadata_pheno$parameter)
  # metadata_pheno_value_size <- length(metadata_pheno$value)
  # if (metadata_pheno_parameter_size != metadata_pheno_value_size) {
  #   stop(paste0("metadata length mismatch: ",
  #               "parameter=", metadata_pheno_parameter_size,
  #               " vs. value=", metadata_pheno_value_size))
  # }
  
  # # --- metadata -> pedigree
  # metadata_pedigree <- data.frame(
  #   parameter = as.vector(c("designation", "mother", "father", "yearOfOrigin")),
  #   value = as.vector(c("germplasmName", "mother", "father", "yearOfOrigin"))
  # )
  
  # --- metadata -> geno
  metadata_geno <- data.frame(
    parameter = as.vector(c("input_format", "ploidity")),
    value = as.vector(c("vcf", 2))
  )
  # # --- modifications -> pheno
  # modifications_pheno <- data.frame(
  #   module = "qaRaw",
  #   analysisId = analysisId,
  #   trait = traits,
  #   reason = "none",
  #   row = NA,
  #   value = NA
  # )
  # row.names(modifications_pheno) <- traits
  
  # --- modifications -> geno
  modifications_geno <- data.frame(reason = c(NA),
                            row = c(NA),
                            col = c(NA),
                            value = c(NA))
  
  modifications_geno$analysisId <- analysisId
  modifications_geno$analysisIdName <- "qa_ebs_mda"
  modifications_geno$module <- "qaGeno"
  
  # --- data -> geno_imp
  imp_freq <- impute_gl(data_geno,
                        ploidity = 2,
                        method = 'frequency')
  
  data_geno_imp <- imp_freq$gl
  
  # --- modifications -> geno_imp
  modifications_geno_imp <- imp_freq$imputation_log$log
  
  
  # --- status
  status <- data.frame(
    module = c("qaRaw","qaGeno"),
    analysisId = rep(analysisId,2),
    analysisIdName = c("qa_ebs_pdm","qa_ebs_mda")
  )
  
  # # --- modeling
  # modeling <- data.frame(
  #   module = "qaRaw",
  #   analysisId = analysisId,
  #   trait = traits,
  #   environment = NA,
  #   parameter = "outlierCoefOutqPheno",
  #   value = NA
  # )
  # row.names(modeling) <- traits
  
  # # --- Create final R object
  result <- list(
    data = list(
      # pheno = data_pheno,  # data.frame
      # pedigree = data_pedigree, # data.frame
      geno = data_geno,  # genlight object
      geno_imp = data_geno_imp  # genlight object
    ),
    metadata = list(
      # pheno = metadata_pheno,  # data.frame
      # pedigree = metadata_pedigree,  # data.frame
      geno = metadata_geno  # data.frame
    ),
    modifications = list(
      # pheno = modifications_pheno  # data.frame
      geno = modifications_geno,  # data.frame
      geno_imp = modifications_geno_imp  # data.frame
    ),
    status = status  # data.frame
    # modeling = modeling  # data.frame
  )
  cat("[RESULT] ")
  cat(str(result))
  # check if folder exist or not
  if (!dir.exists(outputPath)) {
    dir.create(outputPath)
  }
  save(result, file = paste0(outputPath, "/", outputFile, ".RData"))
  return(result)
} # end of staRData fxn
