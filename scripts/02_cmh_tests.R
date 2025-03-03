library(parallel)
library(poolSeq)
library(ACER)
library(qvalue)

ClassicalCmhTest <- function(
    snp_table,
    treatment1,
    gen1 = "01",
    treatment2 = treatment1,
    gen2 = "20") {
  
  # Force evaluation of treatment2
  force(treatment2)
  
  # Changes gen2 if dealing with nBO or nB or B
  if (treatment2 == "nBO" | treatment2 == "nB" | treatment2 == "B") {
    gen2 <- "56"
  }
  
  snp_table_filtered <-
    FilterSamples(
      snp_table = snp_table,
      treatment1 = treatment1,
      gen1 = gen1,
      treatment2 = treatment2,
      gen2 = gen2)
  
  replicates <- CountReplicates(snp_table)
  
  # Picks the columns we will use to build the matrix
  alt1_cols <- paste0("alt_", treatment1, "_rep", replicates, "_gen", gen1)
  alt2_cols <- paste0("alt_", treatment2, "_rep", replicates, "_gen", gen2)
  n1_cols <- paste0("N_", treatment1, "_rep", replicates, "_gen", gen1)
  n2_cols <- paste0("N_", treatment2, "_rep", replicates, "_gen", gen2)
  
  # Coerces the dataset into a matrix of 2x2xreplicates dimensions
  # This uses the parallel package to run the loops in parallel
  p_list <- mclapply(1:nrow(snp_table_filtered), function(line) {
    
    vals <- rbind(
      as.numeric(snp_table_filtered[line, alt1_cols]),
      as.numeric(snp_table_filtered[line, alt2_cols]),
      as.numeric(snp_table_filtered[line, n1_cols]),
      as.numeric(snp_table_filtered[line, n2_cols])
    )
    
    # Check for missing values across all replicates
    if (any(is.na(vals))) return(NA_real_)
    
    # Reshape into a 3D array.
    # Each column of 'vals' corresponds to one replicate:
    # The 2x2 matrix is formed column-wise:
    #   - First column: (v1, v2)
    #   - Second column: (v3, v4)
    matrices <- array(vals, dim = c(2, 2, length(alt1_cols)))
    
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

GetNe <- function(
    snp_table,
    treatment1,
    gen1,
    treatment2,
    gen2,
    t = 20) {
  
  # Ne stands for effective population size. This function uses the methods of ACER to estimate Ne
  # Parameters of the estimateNe function are adjusted for my samples
  alt_cmh1 <- snp_table[grep(paste0("alt_", treatment1, "_rep.._gen", gen1), colnames(snp_table))]
  alt_cmh2 <- snp_table[grep(paste0("alt_", treatment2, "_rep.._gen", gen2), colnames(snp_table))]
  cov_cmh1 <- snp_table[grep(paste0("N_", treatment1, "_rep.._gen", gen1), colnames(snp_table))]
  cov_cmh2 <- snp_table[grep(paste0("N_", treatment2, "_rep.._gen", gen2), colnames(snp_table))]
  
  freq_cmh1 <- alt_cmh1/cov_cmh1
  freq_cmh2 <- alt_cmh2/cov_cmh2
  
  # Creates list of Ne and calculate values for each replicate
  Ne <- c()
  for (i in 1:ncol(freq_cmh1)){
    estimated_Ne <- estimateNe(
      p0 = freq_cmh1[,i],
      pt = freq_cmh2[,i],
      cov0 = cov_cmh1[,i],
      covt = cov_cmh2[,i],
      t = t,
      method = c("P.planII"),
      poolSize = c(100, 100))
    
    # In case I want to save these
    pop_name <- gsub(".*_([A-Z]+_rep\\d+).*", "\\1", colnames(freq_cmh1[i]))
    cat(pop_name, "estimated Ne:", estimated_Ne, "\n")
    Ne <- c(Ne, as.integer(unname(estimated_Ne)))
  }
  
  return(Ne)
}

AdaptedCmhTest <- function(
    snp_table,
    treatment1,
    gen1 = "01",
    treatment2 = treatment1,
    gen2 = "20",
    t = 20,
    Ne) {
  
  # Force evaluation of treatment2
  force(treatment2)
  
  # Changes gen2 if dealing with nBO or nB or B
  if (treatment2 == "nBO" | treatment2 == "nB" | treatment2 == "B") {
    gen2 <- "56"
  }
  
  snp_table_filtered2 <-
    FilterSamples(
      snp_table = snp_table,
      treatment1 = treatment1,
      gen1 = gen1,
      treatment2 = treatment2,
      gen2 = gen2)
  
  replicates <- CountReplicates(snp_table = snp_table_filtered)
  
  freq <- GetFreq(snp_table = snp_table_filtered)
  
  cov <- snp_table_filtered[,grep("^N_", colnames(snp_table_filtered))]
  
  # Checks if both generations are the same. If so, differentiate them
  if (gen1 == gen2) {
    gen1 <- 01
    gen2 <- 02
  }
  
  ## Calculates Ne
  # Ne <-
  #   GetNe(
  #     snp_table = snp_table_filtered,
  #     treatment1 = treatment1,
  #     gen1 = gen1,
  #     treatment2 = treatment2,
  #     gen2 = gen2,
  #     t = t)
  
  pvals <- 
    adapted.cmh.test(
      freq = as.matrix(freq),
      coverage = as.matrix(cov),
      Ne = Ne,
      gen = as.numeric(c(gen1,gen2)),
      repl = 1:length(replicates),
      poolSize = rep(c(100, 100), length(replicates)),
      mincov = 1,
      MeanStart = TRUE,
      IntGen = FALSE,
      TA = FALSE,
      order = 0,
      correct = FALSE,
      RetVal = 0)
  
  return(pvals)
}
