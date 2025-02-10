CalculateFst <- function(snp_table,
                         treatment1,
                         gen1,
                         treatment2,
                         gen2) {
  
  # Based on Evolutionary Genetics by Saetre & Ravinet 
  # Get samples
  samples <- unique(gsub("alt_|N_", "", colnames(snp_table[grep("alt_|N_", colnames(snp_table),)])))
  
  # Calculates replicate count
  replicate_count <- as.numeric(max(unique(gsub("[a-zA-Z]+_rep(..)_gen..", "\\1", samples))))
  
  # Defines frequency dataframe and renames it from "alt_" to "p"
  freq <- snp_table[,grep("^alt_", colnames(snp_table))]/snp_table[,grep("^N_", colnames(snp_table))]
  freq <- freq[,grep(paste("alt_", treatment1, "_rep.._gen", gen1, "|", "alt_", treatment2, "_rep.._gen", gen2,sep = ""), colnames(freq))]
  colnames(freq) <- gsub("alt_", "p_", colnames(freq))
  
  # Initializes freqs dataset to add columns as we calculate them
  freqs <- data.frame(matrix(ncol=0, nrow=nrow(snp_table)))
  
  freqs$p1 <- rowSums(freq[,grep(paste("p_", treatment1, "_rep.._gen", gen1, sep = ""), colnames(freq))])/replicate_count
  freqs$q1 <- 1-freqs$p1
  
  freqs$p2 <- rowSums(freq[,grep(paste("p_", treatment2, "_rep.._gen", gen2, sep = ""), colnames(freq))])/replicate_count
  freqs$q2 <- 1-freqs$p2
  
  freqs$pt <- (freqs$p1 + freqs$p2)/2
  freqs$qt <- (freqs$q1 + freqs$q2)/2
  
  freqs$ht <- 2*freqs$pt*freqs$qt
  
  freqs$hs1 <- 2*freqs$p1*freqs$q1
  freqs$hs2 <- 2*freqs$p2*freqs$q2
  
  freqs$hs <- (freqs$hs1 + freqs$hs2)/2
  
  freqs$fst <- (freqs$ht - freqs$hs)/freqs$ht
  
  return(freqs$fst)
}


