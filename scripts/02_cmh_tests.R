source("scripts/00_useful_functions.R")

library(parallel)
library(poolSeq)
library(ACER)

ClassicalCmhTest <- function(snp_table,
                             treatment1,
                             gen1,
                             treatment2,
                             gen2) {
  
  snp_table_filtered <- FilterSamples(snp_table = snp_table,
                                      treatment1 = treatment1,
                                      gen1 = gen1,
                                      treatment2 = treatment2,
                                      gen2 = gen2)
  
  replicates <- CountReplicates(snp_table = snp_table_filtered,
                                treatment1 = treatment1,
                                gen1 = gen1,
                                treatment2 = treatment2,
                                gen2 = gen2)
  
  # Picks the columns we will use to build the matrix
  alt1_cols <- paste("alt_", treatment1, "_rep", replicates, "_gen", gen1, sep = "")
  alt2_cols <- paste("alt_", treatment2, "_rep", replicates, "_gen", gen2, sep = "")
  n1_cols <- paste("N_", treatment1, "_rep", replicates, "_gen", gen1, sep = "")
  n2_cols <- paste("N_", treatment2, "_rep", replicates, "_gen", gen2, sep = "")
  
  # Coerces the dataset into a matrix of 2x2xreplicates dimensions
  # This uses the parallel package to run the loops in parallel
  p_list <- mclapply(1:nrow(snp_table_filtered), function(line) {
    # Initialize 3D array for contingency tables
    matrices <- array(NA_real_, dim = c(2, 2, length(replicates)))
    
    for (i in seq_along(replicates)) {
      # Extract values and ensure they are numeric
      values <- as.numeric(c(
        snp_table_filtered[line, alt1_cols[i]],
        snp_table_filtered[line, alt2_cols[i]],
        snp_table_filtered[line, n1_cols[i]],
        snp_table_filtered[line, n2_cols[i]]
      ))
      
      # Check for missing values
      if (any(is.na(values))) return(NA_real_)
      
      # Fill the 3D array
      matrices[,,i] <- matrix(values, nrow = 2, ncol = 2)
    }
    
    # Run the Mantel-Haenszel test if the matrix is valid
    test <- tryCatch(mantelhaen.test(matrices), error = function(e) NA_real_)
    test$p.value
  }, mc.cores = max(1, detectCores() - 1))
  
  # Convert results to a numeric vector
  p_list <- unlist(p_list)
  
  # Sometimes we get some NA p-values, which are converted to 1
  # This happens when all samples of one side of the comparison (i.e. all O_gen01 samples) have a frequency of 0
  p_list[is.nan(p_list)] <- 1
  
  return(p_list)
}

getNe <- function(snp_table,
                  treatment1,
                  gen1,
                  treatment2,
                  gen2) {
  
  # Ne stands for effective population size. This function uses the methods of ACER to estimate Ne
  # Parameters of the estimateNe function are adjusted for my samples
  alt_cmh1 <- snp_table[grep(paste("alt_", treatment1, "_rep.._gen", gen1, sep=""), colnames(snp_table))]
  alt_cmh2 <- snp_table[grep(paste("alt_", treatment2, "_rep.._gen", gen2, sep=""), colnames(snp_table))]
  cov_cmh1 <- snp_table[grep(paste("N_", treatment1, "_rep.._gen", gen1, sep=""), colnames(snp_table))]
  cov_cmh2 <- snp_table[grep(paste("N_", treatment2, "_rep.._gen", gen2, sep=""), colnames(snp_table))]
  
  freq_cmh1 <- alt_cmh1/cov_cmh1
  freq_cmh2 <- alt_cmh2/cov_cmh2
  
  # Creates list of Ne and calculate values for each replicate
  Ne <- c()
  for (i in 1:ncol(freq_cmh1)) {
    unlisted_ne <- estimateNe(p0 = freq_cmh1[,i],
                              pt = freq_cmh2[,i],
                              cov0 = cov_cmh1[,i],
                              covt = cov_cmh2[,i],
                              t = 20,
                              method=c("P.planII"),
                              Ncensus=NA,
                              poolSize=c(100, 100))
    
    Ne <- c(Ne, as.integer(unname(unlisted_ne)))
  }
  
  return(Ne)
}

AdaptedCmhTest <- function(snp_table,
                           treatment1,
                           gen1,
                           treatment2,
                           gen2) {
  
  # This function uses the poolSeq package to calculate an CMH test adapted to our genomics use case
  
  snp_table_filtered <- FilterSamples(snp_table = snp_table,
                                      treatment1 = treatment1,
                                      gen1 = gen1,
                                      treatment2 = treatment2,
                                      gen2 = gen2)
  
  replicates <- CountReplicates(snp_table = snp_table_filtered,
                                treatment1 = treatment1,
                                gen1 = gen1,
                                treatment2 = treatment2,
                                gen2 = gen2)
  
  freq <- GetFreq(snp_table = snp_table_filtered,
                  treatment1 = treatment1,
                  gen1 = gen1,
                  treatment2 = treatment2,
                  gen2 = gen2)
  
  cov <- snp_table_filtered[,grep("^N_", colnames(snp_table_filtered))]
  
  
  # Checks if both generations are the same. If so, differentiate them
  if (gen1 == gen2) {
    gen1 <- 01
    gen2 <- 02
  }
  
  # Calculates Ne
  Ne <- getNe(snp_table = snp_table_filtered,
              treatment1 = treatment1,
              gen1 = gen1,
              treatment2 = treatment2,
              gen2 = gen2)
  
  pvals <- adapted.cmh.test(freq = as.matrix(freq),
                            coverage = as.matrix(cov),
                            Ne = Ne,
                            gen = as.numeric(c(gen1,gen2)), # either c(1,20) or c(20,56) - test both 
                            repl = 1:length(replicates),
                            poolSize = rep(100, ncol(freq)),
                            mincov = 1,
                            MeanStart = TRUE,
                            IntGen = FALSE,
                            TA = FALSE,
                            order = 0,
                            correct = FALSE,
                            RetVal = 0)
  
  return(pvals)
}