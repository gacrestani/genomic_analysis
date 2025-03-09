# permutation_tests.R
# Runs the permutation tests using the classic and adapted CMH tests.
#
# inputs: snp_table.rds (x4)
# outputs: perm_pvals.csv - a dataframe with the permutation p-values for each test.

source("scripts/functions.R")

# 0 Initializing ---------------------------------------------------------------
# Read snp tables
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

snp_table_shahrestani_scaled <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani_scaled.rds")

snp_table_regimes_scaled <-
  readRDS("data/processed/processed_snps_abcd_regimes_scaled.rds")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1 Classic --------------------------------------------------------------------
n <- 1000
sapply(
  1:n,
  RunPermutationTestIteration,
  method = "classic",
  filename = "results/permutation_test/classic_OBO01vsOBO20.csv",
  snp_table = snp_table_shahrestani,
  treatment1 = "OBO",
  gen1 = "01",
  treatment2 = "OBO",
  gen2 = "20"
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "classic",
  filename = "results/permutation_test/classic_OB01vsOB20.csv",
  snp_table = snp_table_shahrestani,
  treatment1 = "OB",
  gen1 = "01",
  treatment2 = "OB",
  gen2 = "20"
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "classic",
  filename = "results/permutation_test/classic_nBO01vsnBO56.csv",
  snp_table = snp_table_shahrestani,
  treatment1 = "nBO",
  gen1 = "01",
  treatment2 = "nBO",
  gen2 = "56"
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "classic",
  filename = "results/permutation_test/classic_nB01vsnB56.csv",
  snp_table = snp_table_shahrestani,
  treatment1 = "nB",
  gen1 = "01",
  treatment2 = "nB",
  gen2 = "56"
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "classic",
  filename = "results/permutation_test/classic_O01vsO20.csv",
  snp_table = snp_table_regimes,
  treatment1 = "O",
  gen1 = "01",
  treatment2 = "O",
  gen2 = "20"
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "classic",
  filename = "results/permutation_test/classic_B01vsB56.csv",
  snp_table = snp_table_regimes,
  treatment1 = "B",
  gen1 = "01",
  treatment2 = "B",
  gen2 = "56"
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2 Adapted --------------------------------------------------------------------
n <- 1000
sapply(
  1:n,
  RunPermutationTestIteration,
  method = "adapted",
  filename = "results/permutation_test/adapted_OBO01vsOBO20.csv",
  snp_table = snp_table_shahrestani,
  treatment1 = "OBO",
  gen1 = "01",
  treatment2 = "OBO",
  gen2 = "20",
  Ne = GetNe(
    snp_table = snp_table_shahrestani,
    treatment1 = "OBO",
    gen1 = "01",
    treatment2 = "OBO",
    gen2 = "20")
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "adapted",
  filename = "results/permutation_test/adapted_OB01vsOB20.csv",
  snp_table = snp_table_shahrestani,
  treatment1 = "OB",
  gen1 = "01",
  treatment2 = "OB",
  gen2 = "20",
  Ne = GetNe(
    snp_table = snp_table_shahrestani,
    treatment1 = "OB",
    gen1 = "01",
    treatment2 = "OB",
    gen2 = "20")
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "adapted",
  filename = "results/permutation_test/adapted_nBO01vsnBO56.csv",
  snp_table = snp_table_shahrestani,
  treatment1 = "nBO",
  gen1 = "01",
  treatment2 = "nBO",
  gen2 = "56",
  Ne = GetNe(
    snp_table = snp_table_shahrestani,
    treatment1 = "nBO",
    gen1 = "01",
    treatment2 = "nBO",
    gen2 = "56")
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "adapted",
  filename = "results/permutation_test/adapted_nB01vsnB56.csv",
  snp_table = snp_table_shahrestani,
  treatment1 = "nB",
  gen1 = "01",
  treatment2 = "nB",
  gen2 = "56",
  Ne = GetNe(
    snp_table = snp_table_shahrestani,
    treatment1 = "nBO",
    gen1 = "01",
    treatment2 = "nBO",
    gen2 = "56")
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "adapted",
  filename = "results/permutation_test/adapted_O01vsO20.csv",
  snp_table = snp_table_regimes,
  treatment1 = "O",
  gen1 = "01",
  treatment2 = "O",
  gen2 = "20",
  Ne = GetNe(
    snp_table = snp_table_regimes,
    treatment1 = "O",
    gen1 = "01",
    treatment2 = "O",
    gen2 = "20")
)

sapply(
  1:n,
  RunPermutationTestIteration,
  method = "adapted",
  filename = "results/permutation_test/adapted_B01vsB56.csv",
  snp_table = snp_table_regimes,
  treatment1 = "B",
  gen1 = "01",
  treatment2 = "B",
  gen2 = "56",
  Ne = GetNe(
    snp_table = snp_table_regimes,
    treatment1 = "B",
    gen1 = "01",
    treatment2 = "B",
    gen2 = "56")
)

# ------------------------------------------------------------------------------


# 3 Save results ===============================================================
# Read files and create a new csv with all needed information

obo <- fread("results/permutation_test/classic_OBO01vsOBO20.csv", header = FALSE)
ob  <- fread("results/permutation_test/classic_OB01vsOB20.csv",   header = FALSE)
nbo <- fread("results/permutation_test/classic_nBO01vsnBO56.csv", header = FALSE)
nb  <- fread("results/permutation_test/classic_nB01vsnB56.csv",   header = FALSE)
b   <- fread("results/permutation_test/classic_B01vsB56.csv",     header = FALSE)
o   <- fread("results/permutation_test/classic_O01vsO20.csv",     header = FALSE)

obo_adapted <- fread("results/permutation_test/adapted_OBO01vsOBO20.csv", header = FALSE)
ob_adapted  <- fread("results/permutation_test/adapted_OB01vsOB20.csv",   header = FALSE)
nbo_adapted <- fread("results/permutation_test/adapted_nBO01vsnBO56.csv", header = FALSE)
nb_adapted  <- fread("results/permutation_test/adapted_nB01vsnB56.csv",   header = FALSE)
b_adapted   <- fread("results/permutation_test/adapted_B01vsB56.csv",     header = FALSE)
o_adapted   <- fread("results/permutation_test/adapted_O01vsO20.csv",     header = FALSE)

perm_pvals <-
  cbind(
    obo,
    ob,
    nbo,
    nb,
    o,
    b,
    obo_adapted,
    ob_adapted,
    nbo_adapted,
    nb_adapted,
    o_adapted,
    b_adapted)

colnames(perm_pvals_adapted) <-
  c("obo_adapted",
    "ob_adapted",
    "nbo_adapted",
    "nb_adapted",
    "o_adapted",
    "b_adapted")

perm_pvals <- cbind(perm_pvals, perm_pvals_adapted)

fwrite(perm_pvals, file = "results/perm_pvals.csv")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~