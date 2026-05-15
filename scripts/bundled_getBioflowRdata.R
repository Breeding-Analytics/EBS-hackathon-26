# Auto-generated file. Do not edit directly.
# Source scripts are maintained in modular files under scripts/.
# Generated on: 2026-05-15 08:23:04

# ---- BEGIN: packages_verification.R ----
ensure_cran_packages <- function(packages, repos = "https://cloud.r-project.org") {
  missing <- packages[!vapply(
    packages,
    function(pkg) requireNamespace(pkg, quietly = TRUE),
    FUN.VALUE = logical(1)
  )]

  if (length(missing) > 0) {
    tryCatch(
      install.packages(missing, repos = repos),
      error = function(e) {
        stop(
          sprintf("Failed to install required CRAN packages: %s", paste(missing, collapse = ", ")),
          call. = FALSE
        )
      }
    )
  }

  still_missing <- packages[!vapply(
    packages,
    function(pkg) requireNamespace(pkg, quietly = TRUE),
    FUN.VALUE = logical(1)
  )]

  if (length(still_missing) > 0) {
    stop(
      sprintf("Required CRAN packages are still unavailable: %s", paste(still_missing, collapse = ", ")),
      call. = FALSE
    )
  }
}
# ---- END: packages_verification.R ----

# ---- BEGIN: read_geno_functions.R ----
#' @param geno_metadata Dataframe with the expected columns id, chrom, pos, ref, alt
#'
#' @return Same input dataframe but added a `filter` column if true the marker don't meet
#' the minimal requirements: no duplicated markers, physical position, bi-allelic, 
#'  allelic information, SNP.
#' @export
#'
#' @examples
#' geno_metadata(meta_daya)
process_metadata <- function(geno_metadata){
  # Expected fields id Chrom pos ref alt
  meta_columns <- c('id', 'chrom', 'pos', 'ref', 'alt')
  
  # Check if the columns in the input file match the expected column names
  if (length(intersect(meta_columns, colnames(geno_metadata))) != 5) {
    cli::cli_abort("metadata file doesn't have the expected column names,
                   have: {colnames(geno_metadata)}")
  }
  
  geno_metadata['filter'] <- FALSE
  
  # Flag markers with physical location
  geno_usable_idx <- which(!complete.cases(geno_metadata[,c('chrom','pos')]),)
  #geno_metadata[geno_usable_idx,'filter'] <- TRUE
  no_pos <- length(geno_usable_idx)
  
  if (no_pos > 0){
    cli::cli_inform("Were found {no_pos} markers
                   without physicall location, put into unk chrom and consecutive position")
    
    geno_metadata[geno_usable_idx, 'chrom'] <- "unk"
    geno_metadata[geno_usable_idx, 'pos'] <- seq(1, no_pos)
  }
  
  # Flag colocalized markers
  
  dup_pos_idx <- which(duplicated(geno_metadata[,c('chrom','pos')]))
  no_dup_pos <- length(dup_pos_idx)
  if (no_dup_pos > 0 ){
    cli::cli_inform("Were found {no_dup_pos} markers with duplicated position")
    geno_metadata[dup_pos_idx,'filter'] <- TRUE
  }
  
  # Check duplicated ids
  id_dup <- duplicated(geno_metadata$id)
  
  if(sum(id_dup) > 0){
    cli::cli_inform("Were found {sum(id_dup)} markers with duplicated id. \n
                  Renamed with chrom_pos nomenclature")  
  }
  
  # Count commas in alt column to get allele count
  allele_count <- nchar(geno_metadata$alt) - nchar(gsub(",", "", geno_metadata$alt)) + 1
  
  # Check if ref matches single ACGT pattern
  ref_is_acgt <- grepl("^[ACGT]$", geno_metadata$ref)
  ref_len <- ifelse(ref_is_acgt, nchar(geno_metadata$ref), NA)
  
  # Check if alt matches single ACGT pattern
  alt_is_acgt <- grepl("^[ACGT]$", geno_metadata$alt)
  alt_len <- ifelse(alt_is_acgt, nchar(gsub(",", "", geno_metadata$alt)), NA)
  
  # Update id for duplicates
  row_num <- seq_len(nrow(geno_metadata))
  dup_row_idx <- which(id_dup)
  geno_metadata$id[dup_row_idx] <- paste0(geno_metadata$chrom[dup_row_idx], "_", geno_metadata$pos[dup_row_idx])
  
  # Add new columns
  geno_metadata$allele_count <- allele_count
  geno_metadata$ref_len <- ref_len
  geno_metadata$alt_len <- alt_len
  
  # Flag no reference alleles
  ref_alt_na_idx <- which(is.na(geno_metadata$ref_len) | is.na(geno_metadata$alt_len))
  no_ref_alt <- length(ref_alt_na_idx)
  
  if(no_ref_alt > 0){
    cli::cli_inform("Were found {no_ref_alt} without ref and alt allele data")
    geno_metadata[ref_alt_na_idx,'filter'] <- TRUE
  }
  
  # flag multiallelic 
  multi_allelic_idx <- which(geno_metadata$allele_count > 1)
  multi_allelic <- length(multi_allelic_idx)
  
  if(multi_allelic > 0){
    cli::cli_warn("Were found {multi_allelic} multi-allelic markers")
    geno_metadata[multi_allelic_idx,'filter'] <- TRUE
  }

  # flag indels
  indel_idx <- which(geno_metadata$ref_len > 1)
  indel_idx <- unique(c(indel_idx, which(geno_metadata$alt_len > 1)))
  
  indel <- length(indel_idx)
  
  if(indel > 0){
    cli::cli_inform("Were found {indel} indel markers")
    geno_metadata[indel_idx,'filter'] <- TRUE
  }
  
  return(geno_metadata)
}

#' From locus genotype call data, get the allelic dosage given the alleles and ploidity
#'
#' This function takes a list of genotype calls, a named vector of allele counts,
#' and the ploidity level as input, and returns a list of allelic dosages for the
#' genotype calls. The allelic dosage is the count of the alternative allele in
#' the genotype call.
#'
#' The function first generates all possible genotype calls for the given alleles
#' and ploidity level using the `get_all_gt_calls` function. It then calculates
#' the allelic dosages for these possible genotype calls using the `convert_gt_to_dosage`
#' function, treating the second allele as the alternative allele.
#'
#' Finally, the function replaces the genotype calls in the input list with their
#' corresponding allelic dosages using the `replace_strings_with_integers` function.
#'
#' If a genotype call in the input list is not found in the set of possible genotype
#' calls, its allelic dosage will be set to NA.
#'
#' @param l A list of genotype calls, e.g., c("AG", "GG", "AA").
#' @param alleles ref and alternative allele, e.g., c("A", "C").
#' @param ploidity Integer. The ploidity level of the organism.
#'
#' @return A list of integers, representing the allelic dosages for the input
#'         genotype calls.
#' @export
#'
#' @examples
#' genotypes <- c("AG", "GG", "AA")
#' allele_counts <- c(A = 10, G = 20)
#' get_allelic_dosage(genotypes, allele_counts, 2)  # Returns list(1, 2, 0)
#'
#' # Example with missing genotype call
#' genotypes <- c("AG", "XX", "AA")
#' allele_counts <- c(A = 10, G = 20)
#' get_allelic_dosage(genotypes, allele_counts, 2)  # Returns list(1, NA, 0)
get_allelic_dosage <- function(l, alleles, ploidity, sep = "") {
  alleles_c <- unlist(stringr::str_split(alleles, "/"))
  # All possible genotype calls
  possible_gt_calls <- get_all_gt_calls(alleles_c, ploidity, sep)
  # Get dosage given the alternative allele
  possible_dosage <- convert_gt_to_dosage(possible_gt_calls, alleles_c[2], ploidity,sep)
  dosages <- replace_strings_with_integers(possible_dosage, l)
  return(dosages)
}

#' Get all possible genotype calls given a unique set of alleles
#'
#' This function generates all possible genotype calls for a given set of alleles
#' and ploidity level. The genotype calls are represented as strings of characters,
#' with each allele being a single character.
#'
#' The function uses a recursive approach to generate all possible combinations of
#' alleles for the specified ploidity level. For example, with two alleles "A" and "B",
#' and a ploidity of 2 (diploid), the function would generate the following genotype
#' calls: "AA", "AB", "BA", "BB".
#'
#' @param alleles List[String]. A list of unique alleles, e.g., c("A", "B", "C").
#' @param ploidity Integer. The ploidity level of the organism.
#'
#' @return List[String]. A list of all possible genotype calls for the given alleles
#'         and ploidity level.
#' @export
#'
#' @examples
#' get_all_gt_calls(c("A", "B"), 2)  # Returns c("AA", "AB", "BA", "BB")
#' get_all_gt_calls(c("A", "B", "C"), 3)  # Returns all 27 possible triploid calls
#' get_all_gt_calls(c("A"), 1)  # Returns c("A")
get_all_gt_calls <- function(alleles, ploidity, sep = "") {
  generate_calls <- function(prefix, ploidity, del = sep) {
    if (ploidity == 0) {
      if(nchar(sep) > 0){
        out <- substr(prefix, 1, nchar(prefix)-1)  
      } else {
        out <- substr(prefix, 1, nchar(prefix))  
      }
      
      return(out)
    }
    calls <- c()
    for (allele in alleles) {
      call <- paste0(prefix, allele,del)
      calls <- c(calls, generate_calls(call, ploidity - 1))
    }
    return(calls)
  }
  generate_calls("", ploidity)
}

#' Given a list of genotype calls, get the dosage of each one
#'
#' This function takes a list of genotype calls, an alternative allele, and the ploidity level
#' of the organism as input, and returns a list of allelic dosages corresponding to each
#' genotype call in the input list.
#'
#' @param locus List. A list of genotype calls, e.g., c("AG", "GG", "AA").
#' @param alt_allele String. The alternative allele, e.g., "A", "G".
#' @param ploidity Integer. The ploidity level of the organism, default is 2 (diploid).
#'
#' @return A list of integers, representing the allelic dosages of the alternative allele
#'         for each genotype call in the input list.
#' @export
#'
#' @examples
#' convert_gt_to_dosage(c("AG", "GG", "AA"), "A")  # Returns list(1, 0, 2)
#' convert_gt_to_dosage(c("AG", "GG", "AA"), "G")  # Returns list(1, 2, 0)
#' convert_gt_to_dosage(c("AAA", "GGG"), "A", 3)  # Returns list(3, 0) (triploid)
#' convert_gt_to_dosage(c(NA, "AG"), "A")  # Returns list(NA, 1)
convert_gt_to_dosage <- function(locus, alt_allele, ploidity = 2,sep="") {
  l <- sapply(locus,
              genocall_to_allelic_dosage,
              alt_allele = alt_allele,
              ploidity = ploidity,
              sep=sep)
  return(l)
}


#' Genotype call to allelic dosage of alternative allele
#'
#' This function takes a genotype call, an alternative allele, and the ploidity level
#' of the organism as input, and returns the allelic dosage of the alternative allele
#' in the genotype call.
#'
#' The genotype call is expected to be a string of characters representing the alleles,
#' with each allele being a single character. For example, "AG" represents a diploid
#' genotype with one allele being "A" and the other being "G".
#'
#' The allelic dosage is the count of the alternative allele in the genotype call.
#' For example, if the genotype call is "AG" and the alternative allele is "A", the
#' allelic dosage would be 1.
#'
#' If the genotype call is missing (represented as NA or an empty string), the
#' function returns NA.
#'
#' @param genotype_call String. Genotype call, e.g., "AG", "AAA" (for triploid).
#' @param alt_allele String. Alternative allele, e.g., "A", "G".
#' @param ploidity Integer. Ploidity level of the organism, default is 2 (diploid).
#'
#' @return Integer. The allelic dosage of the given genotype call for the alternative allele.
#' @export
#'
#' @examples
#' genocall_to_allelic_dosage("AG", "A")  # Returns 1
#' genocall_to_allelic_dosage("GG", "A")  # Returns 0
#' genocall_to_allelic_dosage("AAA", "A", 3)  # Returns 3 (triploid)
#' genocall_to_allelic_dosage(NA, "A")  # Returns NA
genocall_to_allelic_dosage <- function(genotype_call, alt_allele, ploidity = 2,sep="") {
  if (!is.na(nchar(genotype_call))) {
    # remove separators (and normalize phasing if present)
    genotype_call <- gsub("\\|", sep, genotype_call)
    if(sep != ""){
      genotype_call <- gsub(sep, "", genotype_call, fixed = TRUE)
    }
    # Genotype call successfully genotyped
    allele_length <- nchar(genotype_call) / ploidity
    
    # List with each allele as element
    split_genotype <- substring(genotype_call, seq(1, nchar(genotype_call), allele_length),
                                seq(allele_length, nchar(genotype_call), allele_length))
    
    # Matches of alt allele are the dosage
    dosage <- length(which(split_genotype == alt_allele))
    return(dosage)
  } else {
    # Genotype call missed
    return(NA)
  }
}

#' Replace a list of strings with their corresponding integer values
#'
#' This function takes two inputs: a named list or vector with string keys and integer values,
#' and a list of strings to be replaced. It replaces each string in the second list with the
#' corresponding integer value from the first list, based on the string keys.
#'
#' If a string in the second list does not have a corresponding key in the first list,
#' it will be replaced with NA.
#'
#' @param lookup_table A named list or vector with string keys and integer values.
#' @param strings_to_replace A list of strings to be replaced with their corresponding integer values.
#'
#' @return A list of integers, where each string in the input list has been replaced with its
#'         corresponding integer value from the lookup table, or NA if no match was found.
#' @export
#'
#' @examples
#' lookup <- c(A = 1, B = 2, C = 3)
#' strings <- c("B", "A", "D", "C")
#' replace_strings_with_integers(lookup, strings)  # Returns list(2, 1, NA, 3)
replace_strings_with_integers <- function(lookup_table, strings_to_replace) {
  # Use match to find the indices of the strings in the lookup table
  indices <- match(strings_to_replace, names(lookup_table))
  
  # Replace the strings with the corresponding integer values
  # or NA if no match was found
  integer_values <- lookup_table[indices]
  
  return(integer_values)
}

get_loc_missing <- function(gl) {
  # Get the number of occurrences of NAs (missing data) for each marker
  NA_counts <- adegenet::glNA(gl)
  
  # Divide the NA counts by the total number of samples to get the missing rate
  NA_counts <- NA_counts / adegenet::nInd(gl)/max(adegenet::ploidy(gl))
  
  return(NA_counts)
}

#' Get Locus Missingness
#'
#' Compute the missing rate for each locus (marker) in a genlight object.
#' The missing rate is a value between 0 and 1, where 0 indicates no missing data
#' for that locus, and 1 indicates that all samples have missing data for that locus.
#'
#' @param gl A genlight object.
#'
#' @return A numeric vector of length equal to the number of loci (markers),
#'   containing the missing rate for each locus.
#'
#' @export
#'
#' @examples
#' data(example_genlight)
#' loc_miss <- get_loc_missing(example_genlight)
#' head(loc_miss)
get_ind_missing <- function(gl) {
  # Convert the genlight object to a matrix and identify missing genotype calls
  mt <- is.na(as.matrix(gl))
  
  # Calculate the proportion of missing data for each individual (row)
  ind_miss <- Matrix::rowSums(mt) / adegenet::nLoc(gl)
  
  return(ind_miss)
}

#' Get Overall Missingness
#'
#' Compute the overall missing rate for a genlight object.
#' The overall missing rate is the proportion of missing genotype calls
#' across all individuals and loci in the dataset.
#'
#' @param gl A genlight object.
#'
#' @return A single numeric value representing the overall missing rate.
#'
#' @export
#'
#' @examples
#' data(example_genlight)
#' overall_miss <- get_overall_missingness(example_genlight)
#' print(overall_miss)
get_overall_missingness <- function(gl) {
  # Convert the genlight object to a matrix
  mt <- as.matrix(gl)
  
  # Identify missing genotype calls
  mt <- is.na(mt)
  
  # Calculate the overall missing rate
  overall_miss <- sum(mt) / (nrow(mt) * ncol(mt))
  
  return(overall_miss)
}

get_heterozygosity_metrics <- function(gl, ploidy = 2){
  # Boolean matrix of genotype calls where 0 > dosage < ploidy
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_loc <- Matrix::colSums(het_ind_loc , na.rm = T)/adegenet::nInd(gl)
  het_ind <- Matrix::rowSums(het_ind_loc , na.rm = T)/adegenet::nLoc(gl)
  return(list(het_ind = het_ind, het_loc = het_loc))
}

get_loc_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_loc <- Matrix::colSums(het_ind_loc , na.rm = T)/adegenet::nInd(gl)
  return(het_loc)
}

get_ind_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_ind <- Matrix::rowSums(het_ind_loc , na.rm = T)/adegenet::nLoc(gl)
  return(het_ind)
}

get_maf <- function(gl, ploidy = 2){
  alf <- adegenet::glMean(gl)*(1/ploidy)
  maf <- ifelse(alf > 0.5, 1 - alf, alf)
  return(maf)
}

get_inbreeding <- function(gl){
  p <- 1 - gl@other$loc.metrics$maf
  he <- 2 * gl@other$loc.metrics$maf * p
  Fis <- 1 - (gl@other$loc.metrics$loc_het/he)
  return(Fis)
}

recalc_metrics <- function(gl){

  gl@other$loc.metrics <- data.frame(
    maf = get_maf(gl, max(adegenet::ploidy(gl))),
    loc_miss = get_loc_missing(gl),
    loc_het = get_loc_heterozygosity(gl, max(adegenet::ploidy(gl)))
  )
  
  # Use already calculated loc stats
  gl@other$loc.metrics$loc_Fis <- get_inbreeding(gl)
  
  gl@other$ind.metrics <- data.frame(
                                     ind_miss = get_ind_missing(gl),
                                     ind_het = get_ind_heterozygosity(gl, max(adegenet::ploidy(gl)))
                                     )
  return(gl)
}


get_overall_summary <- function(gl){
  ninds <- adegenet::nInd(gl)
  nlocs <- adegenet::nLoc(gl)
  overall_missiness <- mean(gl@other$ind.metrics$ind_miss)
  overall_heterozygosity <- mean(gl@other$ind.metrics$ind_het)
  overall_maf <- mean(gl@other$loc.metrics$maf, na.rm = T)
  
  out <- list(
    nind = ninds,
    nloc = nlocs,
    ov_miss = overall_missiness,
    ov_het = overall_heterozygosity,
    ov_maf = overall_maf
  )
  return(out)
}

# Function to infer ploidy from a single genotype call
get_ploidy_from_gt <- function(gt_string) {
  if (is.na(gt_string)) return(NA)
  # Count separators (/ or |) and add 1
  separators <- nchar(gt_string) - nchar(gsub("[/|]", "", gt_string))
  return(separators + 1)
}

#' read_vcf
#'
#' This function reads a VCF file (compressed or uncompressed) and converts it into a genlight object.
#'
#' @param path String. Path to the VCF file. It could be compressed.
#' @param ploidity Integer. Ploidity level of the organism. (Default = 2)
#' @param na_reps Vector. A vector containing the NA representations of genotype calls (default: empty).
#'
#' @return A genlight object.
#' @export
#'
#' @examples
#' fl = "https://github.com/Breeding-Analytics/cgiarGenomics/raw/main/tests/vcf_fmt/diploid.vcf.gz"
#' tempfl <- tempfile(pattern = 'diploid', fileext = '.vcf.gz')
#' download.file(fl, destfile = tempfl)
#' dat.dose.vcf = read_vcf(tempfl, ploidity = 2)
#' print(dat.dose.vcf)
#' plot(dat.dose.vcf)
read_vcf <- function(path, na_reps = c("-", "./."), sep="/") {
  
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  # Read the VCF file
  vcf <- vcfR::read.vcfR(path)
  
  gt_matrix <- vcfR::extract.gt(vcf, return.alleles = TRUE)
  ploidy_matrix <- apply(gt_matrix, c(1, 2), get_ploidy_from_gt)
  # Get the most common ploidy level
  ploidity <- as.numeric(names(sort(table(ploidy_matrix), decreasing = TRUE))[1])

  # Get the metadata from the VCF file
  meta_vcf <- as.data.frame(vcfR::getFIX(vcf))
  
  # Rename columns and select only the needed ones
  meta <- data.frame(
    id = meta_vcf$ID,
    chrom = meta_vcf$CHROM,
    pos = meta_vcf$POS,
    ref = meta_vcf$REF,
    alt = meta_vcf$ALT,
    stringsAsFactors = FALSE
  )
  
  meta <- process_metadata(meta)
  mt <- t(vcfR::extract.gt(vcf, return.alleles = TRUE)[!meta$filter,])
  if (length(na_reps) > 0) {
    idx <- which(mt %in% na_reps)
    mt[idx] <- NA
  }
  
  # check ploidity
  mt_gt_str <- matrix(gsub(sep, "", mt), nrow = dim(mt)[1], ncol = dim(mt)[2])
  gc_len <- apply(mt_gt_str, 2, function(x) max(nchar(x), na.rm = TRUE))
  max_dosage <- max(gc_len, na.rm = TRUE)
  
  if(max_dosage < ploidity){
    cli::cli_warn("Max dosage ({max_dosage}) lower than ploidy lvl ({ploidity})")
  }
  
  if (max_dosage > ploidity){
    cli::cli_abort("Max dosage ({max_dosage}) higher than ploidy lvl ({ploidity})")
  }
  
  individuals <- rownames(mt)
  
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep=sep)

  gt <- mapply(function(col, arg, ploidity, sep) get_allelic_dosage(mt[,col], arg, ploidity, sep),
               col = seq(1, dim(mt)[2]), 
               arg = allele_set,
               ploidity = ploidity,
               sep = sep)
  
  gl <- new("genlight",
            gt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl <- recalc_metrics(gl)
  return(gl)
}

#' Filter function
#' 
#' This function generalizes the filtering functions using the parameter name
#' and comparing in locus or individuals using the provided comparision operator.
#' A list indicating the used thresold, if the filter was performed over individuals
#' or locus and the indices of elements that meet the comparision.
#'
#' @param gl 
#' @param parameter 
#' @param threshold 
#' @param comparison_operator 
#'
#' @return
#' @export
#'
#' @examples
filter_gl <- function(gl, parameter, threshold, comparison_operator){
  
  # Verify if parameter exist on the gl and get the margin (ind, loc)
  filter_margin <- get_parameter_margin(gl, parameter)
  comparison_operator <- match.arg(comparison_operator, choices = c(">", ">=", "<", "<="))
  comparison_func <- match.fun(comparison_operator)
  # Verify if threshold is a float value
  if(threshold > 1 ){
    cli::cli_abort("`threshold`: {threshold} is greather of equal to 1, correct it")
  }
  
  if(filter_margin == 'loc'){
    index <- which(comparison_func(gl@other$loc.metrics[parameter], threshold))
    filter_out <- gl@loc.names[-c(index)]
    
  } else {
    index <- which(comparison_func(gl@other$ind.metrics[parameter], threshold))
    filter_out <- gl@ind.names[-c(index)]
    
  }
  
  out <- list(param = parameter,
              operator = comparison_operator,
              threshold = threshold,
              filter_margin = filter_margin,
              index = index,
              filter_out = filter_out)
  
  return(out)
}

get_parameter_margin <- function(gl, param_name){
  # Get the expected parameters
  loc_metric_names <- colnames(gl@other$loc.metrics)
  ind_metric_names <- colnames(gl@other$ind.metrics)
  
  param_name = match.arg(param_name, choices = c(loc_metric_names, ind_metric_names))
  
  if(param_name %in% loc_metric_names){
    by = 'loc'
  } else {
    by = 'ind'
  }
  return(by)
}

#' Apply a sequence of filterings over a gl object
#'
#' filt_sequence named list (param = param_name, threshold: t, operator: op)
#' @param gl 
#' @param filt_sequence 
#'
#' @return
#' @export
#'
#' @examples
apply_sequence_filtering <- function(gl, filt_sequence){
  if(!rlang::is_bare_list(filt_sequence)){
    cli::cli_abort("Provide a list of filter operations in `filt_sequence`")
  }
  if(!inherits(gl, "genlight")){
    cli::cli_abort("`gl` is not a genlight class")
  }
  
  # Allways add at the end a filter step to remove all NA loci and ind
  
  locNA_filt <- list("loc_miss", "<", 1)
  indNA_filt <- list("ind_miss", "<", 1)
  
  
  allNA_steps <- list(
    locNA_filt,
    indNA_filt
  )
  
  
  filt_NA_seq <- lapply(allNA_steps, function(x){
    setNames(as.list(x), c("param", "operator", "threshold"))
  })
  filt_sequence <- append(filt_sequence, filt_NA_seq)
  # Duplicate the gl object to perform the filtering
  working_gl <- gl
  filtering_log <- list()
  previous_margin <- ""
  
  for (i_step in 1:length(filt_sequence)){
    filt_step <- filt_sequence[[i_step]]
    param <- filt_step[['param']]
    threshold <- filt_step[['threshold']]
    operator <- filt_step[['operator']]
  
    i_filt_out <- filter_gl(working_gl,
                            parameter = param,
                            threshold = threshold,
                            comparison_operator = operator)
    
    i_margin <- get_parameter_margin(gl, param)
    
    if(length(i_filt_out$index) > 0){
      if(i_margin == "loc"){
        working_gl <- working_gl[,i_filt_out$index]
      } else {
        working_gl <- working_gl[i_filt_out$index,]
      }
    }
    working_gl <- recalc_metrics(working_gl)
    filtering_log[[glue::glue("{param}_{i_step}")]] <- i_filt_out
  }
  
  return(list(gl = working_gl, filt_log = filtering_log))
}

get_filter_log <- function(filter_step_log, geno_data){
  print("filter_processing...")
  base_loc_names <- adegenet::locNames(geno_data)
  base_ind_names <- adegenet::indNames(geno_data)
  out <- purrr::map_df(filter_step_log, function(filter_step){
    
    if(length(filter_step$filter_out) > 0){
      
      reason <- paste(filter_step$filter_margin,
                      filter_step$param,
                      filter_step$operator,
                      filter_step$threshold,
                      sep = '_')
      
      if(filter_step$filter_margin == 'loc'){
        loc_idx <- which(filter_step$filter_out %in% base_loc_names)
        col_data <- loc_idx
        row_data <- rep(NA, length(loc_idx))
        
      } else {
        ind_idx <- which(filter_step$filter_out %in% base_ind_names)
        col_data <- rep(NA, length(ind_idx))
        row_data <- ind_idx
      }
      
      filt_step_log <- data.frame(
        reason = rep(reason, length(filter_step$filter_out)),
        row = row_data,
        col = col_data,
        value = rep(NA, length(filter_step$filter_out))
      )
      return(filt_step_log)
    }
  })
  return(out)
}

#' Imputation with allele frequency
#' 
#' Assuming a bi-allelic marker, using the observed allelic frequency for one allele
#' is sampled the genotype call for any ploidity level
#'
#' @param q_frq 
#' @param ploidity 
#'
#' @return
#' @export
#'
#' @examples
i_freq_impute <- function(q_frq, ploidity = 2){
  if(!rlang::is_integerish(ploidity)){
    cli::cli_abort("`ploidity` is not an integer: {ploidity}")  
  }
  dosage <- 0
  for(i_chromatid in seq(ploidity)){
    i_dosage <- sample(c(1,0), size = 1, 
                       prob = c(q_frq, 1 - q_frq), replace = T)
    dosage <- dosage + i_dosage
  }
  return(dosage)
}

freq_impute <- function(gl, mt, ploidity){
  mt <- as.matrix(gl)
  # Get the allelic frequencies
  q_allele <- adegenet::glMean(gl)
  # Linear index of nas
  idx_na <- which(is.na(mt))
  na_loc_idx <- sapply(idx_na, function(x){
    loc_idx <- ceiling(x/nrow(mt))
    return(loc_idx)
  })
  
  if(length(na_loc_idx) > 0){
    imp <- unname(unlist(lapply(q_allele[na_loc_idx],
                                function(x) {
                                  return(as.numeric(i_freq_impute(q_frq = x, ploidity)))})))
    return(split(idx_na, imp))
  } else {
    return(NULL)
  }
}

apply_imputation <- function(mt, imp_dict){
  if (!is.null(imp_dict)){
    for (dosage in names(imp_dict)) {
      # Convert the list name to a numeric value
      num_dosage <- as.numeric(dosage)
      
      # Get the linear indices associated with this value
      idx <- imp_dict[[dosage]]
      
      # Assign the value to these positions in the matrix
      mt[idx] <- num_dosage
    }
  }
  return(mt)
}
#' Impute a gl object
#'  
#' Impute a genlight object using frequency or random forest method.
#' Returns a list with imputed genlight object and imputation log.
#'
#' @param gl genlight object
#' @param ploidity ploidy level (default=2)
#' @param method Imputation method: 'frequency' or 'random_forest' (default='frequency')
#' @param nflank Number of flanking markers for RF method (default=100)
#' @param ntree Number of trees for RF method (default=100)
#' @param seed Optional seed for reproducibility (RF method)
#'
#' @return List with elements:
#'   - gl: imputed genlight object
#'   - log: imputation dictionary (imputed positions grouped by dosage)
#' @export
#'
#' @examples
impute_gl <- function(gl, ploidity = 2, method = 'frequency', nflank = 100, ntree = 100, seed = NULL){
  
  loci_all_nas <- adegenet::glNA(gl)/ploidity == adegenet::nInd(gl)
  
  if(sum(loci_all_nas) > 0){
    cli::cli_warn("There are {sum(loci_all_nas)} loci with all missing data")
    # Filter out the all na loci
    all_notna_idxs <- which(!loci_all_nas)
    gl <- gl[,all_notna_idxs]
    mt <- as.matrix(gl)
  }
  
  
  nas_number <- sum(adegenet::glNA(gl))/ploidity
  number_imputations <- nas_number - (sum(loci_all_nas) * adegenet::nInd(gl))
  
  mt <- as.matrix(gl)
  
  
  cli::cli_inform("Missing genotype calls {number_imputations}")
  
  if(method == 'frequency'){
    imp_dict <- freq_impute(gl, mt, ploidity)
  } else if(method == 'random_forest'){
    cli::cli_inform("Imputing with Random Forest (nflank={nflank}, ntree={ntree})")
    imp_dict <- rf_impute(gl, nflank = nflank, ntree = ntree, seed = seed)
  } else {
    cli::cli_abort("Unknown imputation method: {method}. Use 'frequency' or 'random_forest'")
  }
  
  # apply the imputation creating a new gl instance
  imp_mt <- apply_imputation(mt, imp_dict)
  
  
  imp_gl <- new("genlight",
            imp_mt,
            ploidy = ploidity,
            loc.names = gl@loc.names,
            ind.names = gl@ind.names,
            chromosome = gl@chromosome,
            position = gl@position)
  
  adegenet::alleles(imp_gl) <- adegenet::alleles(gl)
  imp_gl <- recalc_metrics(imp_gl)
  
  return(list(gl = imp_gl, log = imp_dict))
}
# ---- END: read_geno_functions.R ----

# ---- BEGIN: getBioflowRdata.R ----
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

getBioflowRData <- function(phenotypeFile, pedigreeFile = NULL, genotypeFile = NULL, 
                     traits, outputPath, outputFile, requestId = NULL) {
  
  # Now we can use ensure_cran_packages since packages_verification.R has been sourced
  ensure_cran_packages(c("vcfR", "adegenet", "cli", "rlang"))
  analysisId <- round(as.numeric(Sys.time()), 0)
  cat(paste0("analysisId: ", analysisId, "\n"))
  
  traits <- unlist(traits)  # ensure object is a vector
  cat(paste0("trait: ", traits, "\n"))
  cat("\n")
  
  # --- data -> pheno
  data_pheno <- read.csv(phenotypeFile, encoding = 'utf-8')

  # # remove non-alphanumeric character in occurrenceName and revise how enviroment is created
  data_pheno$environment <- paste0("Env", paste(data_pheno$year, data_pheno$season, gsub("[^a-zA-Z0-9]", "", data_pheno$occurrenceName), sep = "_"))
  
  # check if variable are present
  if (any(is.na(data_pheno$design)) || any(data_pheno$design == "")) {
    stop("Design information are missing for some or all rows.")
  }
  
  # check design parameter base, rep in EBS pertains to the number of times that entry appear in the occurrence 
  data_pheno[data_pheno$design == "Partially Replicated", "rep"] <- NA
  data_pheno[data_pheno$design == "Augmented", "rep"] <- NA
  data_pheno[data_pheno$design == "Augmented RCBD", "rep"] <- NA # assumes that blockNumber is not NA
  
  # --- data -> pedigree
  data_pedigree <- data.frame(
    germplasmName = as.vector(unique(data_pheno$germplasmName)),
    mother = NA,
    father = NA,
    yearOfOrigin = NA
  )
  
  # --- data -> geno
  data_geno <- read_vcf(genotypeFile) # todo generalize the ploidity level
  # --- metadata -> pheno
  # NOTE: rep() is bugged with python
  # metadata_pheno <- data.frame(
  #   parameter = as.vector(c("stage", "year", "season", "location", "trial", "study", "rep", "iBlock", "row", "col", "designation", "gid", "entryType", rep("trait", nTrait), "environment", "source")),
  #   value = as.vector(c("breedingStage", "year", "season", "site", "experiment", "occurrenceName", "rep", "blockNumber", "paY", "paX", "germplasmName", "germplasmDbId", "entryType", traits, "environment", "ebs-ba"))
  # )
  
  metadata_pheno1 <- data.frame(
    parameter = c("stage", "year", "season", "location", "trial", "study", "rep", "iBlock", "row", "col", "designation", "gid", "entryType"),
    value = c("breedingStage", "year", "season", "site", "experiment", "occurrenceName", "rep", "blockNumber", "paY", "paX", "germplasmName", "germplasmDbId", "entryType")
  )
  metadata_pheno2 <- data.frame(
    parameter = "trait",
    value = traits
  )
  if (is.null(requestId)) {
    metadata_pheno3 <- data.frame(
      parameter = c("environment", "source", "sourceId"),
      value = c("environment", "ebs-ba", NA)
    )
  } else {
    metadata_pheno3 <- data.frame(
      parameter = c("environment", "source", "sourceId"),
      value = c("environment", "ebs-ba", requestId)
    )  
  }
  
  metadata_pheno <- rbind(rbind(metadata_pheno1, metadata_pheno2), metadata_pheno3)
  
  metadata_pheno_parameter_size <- length(metadata_pheno$parameter)
  metadata_pheno_value_size <- length(metadata_pheno$value)
  if (metadata_pheno_parameter_size != metadata_pheno_value_size) {
    stop(paste0("metadata length mismatch: ",
                "parameter=", metadata_pheno_parameter_size,
                " vs. value=", metadata_pheno_value_size))
  }
  
  # --- metadata -> pedigree
  metadata_pedigree <- data.frame(
    parameter = as.vector(c("designation", "mother", "father", "yearOfOrigin")),
    value = as.vector(c("germplasmName", "mother", "father", "yearOfOrigin"))
  )
  
  # --- metadata -> geno
  metadata_geno <- data.frame(
    parameter = as.vector(c("input_format", "ploidity")),
    value = as.vector(c("vcf", 2))
  )
  # --- modifications -> pheno
  modifications_pheno <- data.frame(
    module = "qaRaw",
    analysisId = analysisId,
    trait = traits,
    reason = "none",
    row = NA,
    value = NA
  )
  row.names(modifications_pheno) <- traits
  
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
    analysisId = c(analysisId,analysisId),
    analysisIdName = c("qa_ebs_pdm","qa_ebs_mda")
  )
  
  # --- modeling
  modeling <- data.frame(
    module = "qaRaw",
    analysisId = analysisId,
    trait = traits,
    environment = NA,
    parameter = "outlierCoefOutqPheno",
    value = NA
  )
  row.names(modeling) <- traits
  
  # --- Create final R object
  result <- list(
    data = list(
      pheno = data_pheno,  # data.frame
      pedigree = data_pedigree, # data.frame
      geno = data_geno,  # genlight object
      geno_imp = data_geno_imp  # genlight object
    ),
    metadata = list(
      pheno = metadata_pheno,  # data.frame
      pedigree = metadata_pedigree,  # data.frame
      geno = metadata_geno  # data.frame
    ),
    modifications = list(
      pheno = modifications_pheno, # data.frame
      geno = modifications_geno,  # data.frame
      geno_imp = modifications_geno_imp  # data.frame
    ),
    status = status,  # data.frame
    modeling = modeling  # data.frame
  )
  cat("[RESULT] ")
  cat(str(result))
  # check if folder exist or not
  if (!dir.exists(outputPath)) {
    dir.create(outputPath)
  }
  save(result, file = paste0(outputPath, "/", outputFile, ".RData"))
} # end of staRData fxn
# ---- END: getBioflowRdata.R ----

