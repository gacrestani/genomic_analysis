source("scripts/00_useful_functions.R")

RunPermutationTestIteration <- 
  function(method = c("classic", "adapted"),
           filename,
           iter,
           snp_table,
           treatment1,
           gen1,
           treatment2,
           gen2) {
  
  # WARNING!
  # I think you don't have to run the permutation test if using the adapted 
  # method. Double check literature!
  
  # Assures a correct method is selected
  method <- match.arg(method)
  
  # Measure time needed to run one iteration
  start_time <- Sys.time()
  
  # Filter df to contain only interesting columns
  snp_table_filtered <-
    FilterSamples(snp_table,
                  treatment1 = treatment1,
                  gen1 = gen1,
                  treatment2 = treatment2,
                  gen2 = gen2)
  
  # Define possible choices for the sampling process
  choices_pattern_a <- paste("alt_", treatment1, "_rep.._gen", gen1, sep="")
  choices_pattern_b <- paste("alt_", treatment2, "_rep.._gen", gen2, sep="")
  
  choices_a <- 
    colnames(
      snp_table_filtered[,grep(choices_pattern_a, colnames(snp_table_filtered))]
      )
  
  choices_b <- colnames(
    snp_table_filtered[,grep(choices_pattern_b, colnames(snp_table_filtered))]
    )
  
  choices <- c(choices_a, choices_b)
  
  # Sample (randomize) choices
  choices_rand_alt <- sample(choices, length(choices))
  choices_rand_N <- gsub("alt_", "N_", choices_rand_alt)
  
  choices_rand <- c(rbind(choices_rand_alt, choices_rand_N))
  
  colnames(snp_table_filtered) <- choices_rand
  
  if (method == "classic") {
    
    pval_list <- 
      ClassicalCmhTest(
        snp_table = snp_table,
        treatment1 = treatment1,
        gen1 = gen1,
        treatment2 = treatment2,
        gen2 = gen2)
    
  } else if (method == "adapted") {
    
    pval_list <- 
      AdaptedCmhTest(
        snp_table = snp_table_filtered,
        treatment1 = treatment1,
        gen1 = gen1,
        treatment2 = treatment2,
        gen2 = gen2)
    
  }
  
  pval_list_log <- -log10(pval_list)
  max_pval <- max(pval_list_log)
  
  end_time <- Sys.time()
  end_time-start_time
  print(paste("Iter", iter, "- Time elapsed:", end_time-start_time, "seconds.", "p-val:", max_pval))
  
  fwrite(x = as.list(max_pval),
         file = filename,
         append=TRUE)
  
  return(max_pval)
}