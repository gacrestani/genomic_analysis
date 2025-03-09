# data_preparation.R
# This script prepares the data for the downstream analysis.
# It adds ABS_POS, filters MAF, filters minimum and maximum coverage, and filters chromosomes.
#
# inputs: raw filtered_snps.txt files coming out of genomics pipeline
# outputs: snp_table.Rds - two pairs of SNP tables, one for each naming convention. 
#   The firs pair is untransformed data, the second pair is scaled data.

source("scripts/functions.R")

# 1 Data cleaning --------------------------------------------------------------
snp_table_shahrestani <- as.data.frame(ReadAndPrepare(mode = "shahrestani"))
saveRDS(snp_table_shahrestani,
        "data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <- as.data.frame(ReadAndPrepare(mode = "regimes"))
saveRDS(snp_table_regimes,
        "data/processed/processed_snps_abcd_regimes.rds")

# Create scaled snp tables
snp_table_shahrestani_scaled <- ScaleSnpTable(snp_table_shahrestani)
saveRDS(snp_table_shahrestani_scaled,
        "data/processed/processed_snps_abcd_shahrestani_scaled.rds")

snp_table_regimes_scaled <- ScaleSnpTable(snp_table_regimes)
saveRDS(snp_table_regimes_scaled,
        "data/processed/processed_snps_abcd_regimes_scaled.rds")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
