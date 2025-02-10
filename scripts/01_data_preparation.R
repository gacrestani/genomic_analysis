library(dplyr)
library(data.table)

ReadSnpTable <- function(path = "data/snp_tables/filtered_snps_abcd.txt",
                         mode = "shahrestani") {
  
  snp_table <- data.table::fread(file = path, header = TRUE)
  snp_table <- as.data.frame(snp_table)
  
  # Replaces "." with 0 and converts columns 6 to the end (i.e. those with data) to numeric.
  # TO DO: why are we getting those dots?
  snp_table[snp_table == "."] <- 0
  snp_table[, 6:ncol(snp_table)] <- lapply(snp_table[, 6:ncol(snp_table)], as.numeric)
  
  # Determines the sample naming convention used in my project
  # default:     EBO, EB, CBO, CB
  # shahrestani: OBO, OB, nBO, nB
  # regimes:     O_1~10, B_1~10
  
  if (mode == "regimes") {
    colnames <- colnames(snp_table)
    colnames <- gsub("CBO", "B", colnames(snp_table))
    colnames <- gsub("CB_rep01", "B_rep06", colnames)
    colnames <- gsub("CB_rep02", "B_rep07", colnames)
    colnames <- gsub("CB_rep03", "B_rep08", colnames)
    colnames <- gsub("CB_rep04", "B_rep09", colnames)
    colnames <- gsub("CB_rep05", "B_rep10", colnames)
    
    colnames <- gsub("EBO", "O", colnames)
    colnames <- gsub("EB_rep01", "O_rep06", colnames)
    colnames <- gsub("EB_rep02", "O_rep07", colnames)
    colnames <- gsub("EB_rep03", "O_rep08", colnames)
    colnames <- gsub("EB_rep04", "O_rep09", colnames)
    colnames <- gsub("EB_rep05", "O_rep10", colnames)
    
    colnames(snp_table) <- colnames
    
  } else if (mode == "shahrestani") {
    colnames <- colnames(snp_table)
    
    colnames <- gsub("EB", "OB", colnames)
    colnames <- gsub("CB", "nB", colnames)
    
    colnames(snp_table) <- colnames
  }
  
  return(snp_table)
}

AddAbsPosToSnpTable <- function(snp_table) {
  
  # Drosophila chromosome size data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4352887/ (table 1)
  chromosome_sizes <- c("2L" = 23513712,
                        "2R" = 25286936,
                        "3L" = 28110227,
                        "3R" = 32079331,
                        "X" = 23542271) 
  
  snp_table <- snp_table %>%
    mutate(ABS_POS = case_when(CHROM == "2L" ~ POS,
                               CHROM == "2R" ~ POS + sum(chromosome_sizes[1]),
                               CHROM == "3L" ~ POS + sum(chromosome_sizes[1:2]),
                               CHROM == "3R" ~ POS + sum(chromosome_sizes[1:3]),
                               CHROM == "X" ~ POS  + sum(chromosome_sizes[1:4])))
  
  return(snp_table)
}

FilterMAF <- function(snp_table,
                      limit = 0.001) {
  
  # MAF stands for minor allele frequency
  # If all samples are below a MAF cutoff, the SNP is removed
  samples <- gsub("N_", "", colnames(snp_table)[grepl("^N_", colnames(snp_table))])
  alt_sums <- rowSums(snp_table[paste("alt_", samples, sep="")])
  total_sums <- rowSums(snp_table[paste("N_", samples, sep="")])
  freq_rowsums <- alt_sums / total_sums
  snp_table <- snp_table[freq_rowsums >= limit & freq_rowsums <= (1 - limit),]
  
  return(snp_table)
}

FilterMinAndMaxCov <- function(snp_table,
                               min_cov = 30,
                               max_cov = 500) {
  
  snp_table <- snp_table %>%
    filter(if_all(starts_with("N_"), ~ . >= min_cov & . <= max_cov))
  
  return(snp_table)
}

FilterChromosomes <- function(snp_table = snp_table,
                              chromosomes = c("2L", "2R", "3L", "3R", "X")) {

  snp_table <- snp_table[snp_table$CHROM %in% chromosomes, ]
  
  return(snp_table)
}

ReadAndPrepare <- function(path = "data/snp_tables/filtered_snps_abcd.txt",
                           mode = "shahrestani",
                           limit = 0.001,
                           min_cov = 30,
                           max_cov = 500,
                           chromosomes = c("2L", "2R", "3L", "3R", "X")) {
  
  snp_table <- ReadSnpTable(path, mode)
  snp_table <- AddAbsPosToSnpTable(snp_table)
  snp_table <- FilterMAF(snp_table, limit)
  snp_table <- FilterMinAndMaxCov(snp_table, min_cov, max_cov)
  snp_table <- FilterChromosomes(snp_table, chromosomes)
  
  return(snp_table)
  
}