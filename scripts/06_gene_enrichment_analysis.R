library(biomaRt)

FilterOutFixedSnps <- function(snp_table) {
  
  o_type_samples <- "alt_OB_rep0._gen01|alt_OBO_rep0._gen01"
  
  selected_columns <-
    colnames(snp_table)[grep(o_type_samples, colnames(snp_table))]
  
  rows_to_exclude <-
    apply(freq[selected_columns],
          1,
          function(row) {any(row == 0) || any(row == 1)})
  
  filtered_snp_table <- snp_table[!rows_to_exclude, ]
  
  rows_lost <- sum(rows_to_exclude)
  
  message <-
    paste("Number of SNPs discarded:", format(rows_lost, big.mark = ","))
  
  print(message)
  
  return(filtered_snp_table)
}

GetLinkageWindows <- function(pos, window_size = 50000) {
  return(c(pos - window_size, pos + window_size))
}

DiagnoseSnps <- function(filtered_snp_table) {
  for (i in 1:nrow(filtered_snp_table)) {
    
    GetAllelicTrajectoryPlot(snp_table = filtered_snp_table,
                             snp = paste(filtered_snp_table$CHROM[i], ":", filtered_snp_table$POS[i], sep = ""),
                             treatments = c("OBO", "OB"))
    
    cat("Row ", i, "\n")
    
    # Prompt user input
    response <- readline(prompt = "Press enter to continue. Type anything to exit.\n")
    
    # Check response
    if (toupper(response) != "") {
      cat("Exiting.\n")
      break
    }
    
    cat("Continuing to the next SNP.\n")
  }
}







#gene_list_backup <- gene_list
# gene_list <- gene_list_backup
# 
# filtered_vector <- my_vector[!grepl(pattern, my_vector)]
# 
# gene_list <- gene_list[!grepl("^lncRNA*", gene_list)]
# gene_list_unique <- unique(gene_list)
# 
# cat(gene_list_unique, sep = "\n")

GetGeneList <-
  function(GO_dataframe,
           peaks,
           ensembl,
           window = 10000) {
  
  gene_list <- c()
  
  for (i in 1:length(peaks)) {
    chromosome <- GO_dataframe$CHROM[GO_dataframe$coordinate == peaks[i]]
    pos <- GO_dataframe$POS[GO_dataframe$coordinate == peaks[i]]
    
    result <- getBM(
      attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
      filters = c("chromosome_name", "start", "end"),
      values = list(chromosome, pos-window, pos+window),
      mart = ensembl
    )
    
    print(result)
    gene_list <- c(gene_list, result$external_gene_name)
  }

  return(gene_list)
  
}