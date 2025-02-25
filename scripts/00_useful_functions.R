FilterSamples <- function(snp_table,
                          treatment1,
                          gen1,
                          treatment2,
                          gen2) {
  
  sample1 <- paste(treatment1, "_rep.._gen", gen1, sep = "")
  sample2 <- paste(treatment2, "_rep.._gen", gen2, sep = "")
  
  snp_table_filtered <- snp_table[grep(paste(sample1, "|", sample2, sep = ""),
                                       colnames(snp_table))]
  
  return(snp_table_filtered)
}

CountReplicates <- function(snp_table) {
  
  freq <- GetFreq(snp_table)
  replicates <- unique(gsub("[a-zA-Z]+_[a-zA-Z]+_rep(..)_gen..",
                            "\\1",
                            colnames(freq)))
  
  return(replicates)
}

GetFreq <- function(snp_table) {
  
  alt <- snp_table[,grep("^alt_", colnames(snp_table))]
  cov <- snp_table[,grep("^N_", colnames(snp_table))]
  
  freq <- alt / cov
  
  return(freq)
}
