source("scripts/functions.R")

RunPermutationTestIteration <- function(iter,
                                        snp_table,
                                        treatment1,
                                        gen1,
                                        treatment2,
                                        gen2) {
  
  # Measure time
  start_time <- Sys.time()
  
  # Filter df to contain only interesting columns
  df_filtered <- snp_table[grep(paste(treatment1, "_rep.._gen", gen1, "|", treatment2, "_rep.._gen", gen2, sep = ""), colnames(snp_table))]
  
  # Define possible choices for the sampling process
  choices_pattern_a <- paste("alt_", treatment1, "_rep.._gen", gen1, sep="")
  choices_pattern_b <- paste("alt_", treatment2, "_rep.._gen", gen2, sep="")
  
  choices_a <- colnames(df_filtered[,grep(choices_pattern_a, colnames(df_filtered))])
  choices_b <- colnames(df_filtered[,grep(choices_pattern_b, colnames(df_filtered))])
  
  choices <- c(choices_a, choices_b)
  
  # Sample (randomize) choices
  choices_rand_alt <- sample(choices, length(choices))
  choices_rand_N <- gsub("alt_", "N_", choices_rand_alt)
  
  choices_rand <- c(rbind(choices_rand_alt, choices_rand_N))
  
  colnames(df_filtered) <- choices_rand
  
  pval_list <- ClassicalCmhTest(snp_table = snp_table,
                                 treatment1 = treatment1,
                                 gen1 = gen1,
                                 treatment2 = treatment2,
                                 gen2 = gen2)
  
  pval_list_log <- -log10(pval_list)
  max_pval <- max(pval_list_log)
  
  end_time <- Sys.time()
  end_time-start_time
  print(paste("Iter", iter, "- Time elapsed:", end_time-start_time, "seconds.", "p-val:", max_pval))
  
  fwrite(x = as.list(max_pval),
         file = paste("results/permutation_test/", treatment1, gen1, "vs", treatment2, gen2, ".csv", sep=""),
         append=TRUE)
  
  return(max_pval)
}

snp_table_sharestani <- fread("data/processed/processed_snps_abcd_shahrestani.csv")
snp_table_regimes <- fread("data/processed/processed_snps_abcd_regimes.csv")

n <- 1000
sapply(1:n,
       runPermutationTestIteration,
       snp_table = snptable,
       treatment1 = "OBO",
       gen1 = "01",
       treatment2 = "OBO",
       gen2 = "20")

sapply(1:n,
       runPermutationTestIteration,
       snp_table=snptable,
       treatment1="OB",
       gen1="01",
       treatment2="OB",
       gen2="20")

sapply(1:n,
       runPermutationTestIteration,
       snp_table=snptable,
       treatment1="nBO",
       gen1="01",
       treatment2="nBO",
       gen2="56")

sapply(1:n,
       runPermutationTestIteration,
       snp_table=snptable,
       treatment1="nB",
       gen1="01",
       treatment2="nB",
       gen2="56")



sapply(1:n,
       runPermutationTestIteration,
       snp_table=snptable,
       treatment1="O",
       gen1="01",
       treatment2="O",
       gen2="20")

sapply(1:n,
       runPermutationTestIteration,
       snp_table=snptable,
       treatment1="B",
       gen1="01",
       treatment2="B",
       gen2="56")

# Read files and create a new csv with all needed information
obo <- read.csv("results/OBO01vsOBO20.csv", header = FALSE)
ob <- read.csv("results/OB01vsOB20.csv", header = FALSE)
nbo <- read.csv("results/nBO01vsnBO56.csv", header = FALSE)
nb <- read.csv("results/nB01vsnB56.csv", header = FALSE)
o <- read.csv("results/O01vsO20.csv", header = FALSE)
b <- read.csv("results/B01vsB56.csv", header = FALSE)

perm.pvals <- cbind(obo,ob,nbo,nb,o,b)
colnames(perm.pvals) <- c("obo", "ob", "nbo", "nb", "o", "b")

write.csv(perm.pvals, file = "results/perm_pvals.csv")


perm <- read.csv("results/perm_pvals.csv")

hist(perm$obo)
hist(perm$ob)
hist(perm$nbo)
hist(perm$nb)
hist(perm$o)
hist(perm$b)


## Adapted cmh test

runPermutationTestIteration.adapted <- function(iter,
                                                snp_table,
                                                treatment1,
                                                gen1,
                                                treatment2,
                                                gen2) {
  
  # Measure time
  start.time <- Sys.time()
  #print("Starting...")
  #print(start.time)
  
  # Filter df to contain only interesting columns
  df.filtered <- snp_table[grep(paste(treatment1, "_rep.._gen", gen1, "|", treatment2, "_rep.._gen", gen2, sep = ""), colnames(snp_table))]
  
  # Define possible choices for the sampling process
  choices.pattern.a <- paste("alt_", treatment1, "_rep.._gen", gen1, sep="")
  choices.pattern.b <- paste("alt_", treatment2, "_rep.._gen", gen2, sep="")
  
  choices.a <- colnames(df.filtered[,grep(choices.pattern.a, colnames(df.filtered))])
  choices.b <- colnames(df.filtered[,grep(choices.pattern.b, colnames(df.filtered))])
  
  choices <- c(choices.a, choices.b)
  
  # Sample (randomize) choices
  choices.rand.alt <- sample(choices, length(choices))
  choices.rand.N <- gsub("alt_", "N_", choices.rand.alt)
  
  choices.rand <- c(rbind(choices.rand.alt, choices.rand.N))
  
  colnames(df.filtered) <- choices.rand
  
  p.val.list <- adaptedCMH(snp_table = df.filtered,
                           treatment1 = treatment1,
                           gen1 = gen1,
                           treatment2 = treatment2,
                           gen2 = gen2)
  
  max.p.val <- max(p.val.list)
  end.time <- Sys.time()
  end.time-start.time
  print(paste("Iter", iter, "- Time elapsed:", end.time-start.time, "seconds.", "p-val:", max.p.val))
  
  fwrite(x = as.list(max.p.val),
         file = paste("results/",treatment1,gen1,"vs",treatment2,gen2,"_adapted.csv", sep=""),
         append=TRUE)
  
  return(max.p.val)
}

n <- 1000
sapply(1:n,
       runPermutationTestIteration.adapted,
       snp_table=snptable,
       treatment1="OBO",
       gen1="01",
       treatment2="OBO",
       gen2="20")
