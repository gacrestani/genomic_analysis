# Header =======================================================================
# This script is the main script that runs the entire analysis. It reads the raw
# snp_table, processes it, and runs the CMH tests. The results are saved in the
# results folder.

# 0 - Load useful functions ====================================================
source("scripts/00_useful_functions.R")
source("scripts/01_data_preparation.R")
source("scripts/02_cmh_tests.R")
source("scripts/03_permutation_test.R")
source("scripts/04_plot_functions.R")
source("scripts/05_pca_analysis.R")
source("scripts/06_gene_enrichment_analysis.R")

options(scipen=999) # Disable scientific notation

# 1 - Process raw data =========================================================
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

# 2 - CMH tests ================================================================
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

# Initiates an empty dataframe with the same number of rows as our snp tables
cmh_pvals <-
  data.frame(matrix(NA,
                    nrow = nrow(snp_table_shahrestani),
                    ncol = 0))

cmh_pvals$ABS_POS <- snp_table_shahrestani$ABS_POS
cmh_pvals$CHROM <- snp_table_shahrestani$CHROM

# 2.1 - Unscaled data ==========================================================
# NOTE: These functions take a while to run (even running in parallel)

treatments <- c("OBO", "OB", "nBO", "nB", "O", "B")
gen2 <- c("20", "20", "56", "56", "20", "56")

# Classic CMH
# Takes a while to run
for (i in 1:length(treatments)) {
  
  # Define if using regimes or shahrestani snp table
  ifelse(
    treatments[i] == "O" | treatments[i] == "B",
    snp_table_iter <- snp_table_regimes,
    snp_table_iter <- snp_table_shahrestani
  )
  
  # Run test
  cmh_pvals[[paste0("cmh_classic_", treatments[i])]] <-
    ClassicalCmhTest(
      snp_table = snp_table_iter,
      treatment1 = treatments[i])
  
  # Print test name completed
  print(paste("Classic CMH test done:"), treatment)
  
}

# Adapted CMH
# These are faster than the classical tests
for (i in 1:length(treatments)) {
  
  # Define if using regimes or shahrestani snp table
  ifelse(
    treatments[i] == "O" | treatments[i] == "B",
    snp_table_iter <- snp_table_regimes,
    snp_table_iter <- snp_table_shahrestani
  )
  
  # Run test
  cmh_pvals[[paste0("cmh_adapted_", treatments[i])]] <-
    AdaptedCmhTest(
      snp_table = snp_table_iter,
      treatment1 = treatments[i])
  
  # Print test name completed
  print(paste("Classic CMH test done:", treatments[i]))
  
}


# 2.2 - Scaled data ============================================================
snp_table_shahrestani_scaled <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani_scaled.rds")

snp_table_regimes_scaled <-
  readRDS("data/processed/processed_snps_abcd_regimes_scaled.rds")

treatments <- c("OBO", "OB", "nBO", "nB", "O", "B")
gen2 <- c("20", "20", "56", "56", "20", "56")

# Classic CMH
# Takes a while to run
for (i in 1:length(treatments)) {
  
  # Define if using regimes or shahrestani snp table
  ifelse(
    treatments[i] == "O" | treatments[i] == "B",
    snp_table_iter <- snp_table_regimes_scaled,
    snp_table_iter <- snp_table_shahrestani_scaled
  )
  
  # Run test
  cmh_pvals[[paste0("cmh_classic_scaled_", treatments[i])]] <-
    ClassicalCmhTest(
      snp_table = snp_table_iter,
      treatment1 = treatments[i])
  
  # Print test name completed
  print(paste("Classic CMH test done:"), treatment)
  
}

# Adapted CMH
# These are faster than the classical tests
for (i in 1:length(treatments)) {
  
  # Define if using regimes or shahrestani snp table
  ifelse(
    treatments[i] == "O" | treatments[i] == "B",
    snp_table_iter <- snp_table_regimes,
    snp_table_iter <- snp_table_shahrestani
  )
  
  # Run test
  cmh_pvals[[paste0("cmh_adapted_scaled_", treatments[i])]] <-
    AdaptedCmhTest(
      snp_table = snp_table_iter,
      treatment1 = treatments[i])
  
  # Print test name completed
  print(paste("Classic CMH test done:", treatments[i]))
  
}


# 3 - Permutation tests ========================================================
snp_table_shahrestani <-
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

# 3.1 - Classic ################################################################
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

# Read files and create a new csv with all needed information
obo <- fread("results/permutation_test/classic_OBO01vsOBO20.csv", header = FALSE)
ob  <- fread("results/permutation_test/classic_OB01vsOB20.csv",   header = FALSE)
nbo <- fread("results/permutation_test/classic_nBO01vsnBO56.csv", header = FALSE)
nb  <- fread("results/permutation_test/classic_nB01vsnB56.csv",   header = FALSE)
b   <- fread("results/permutation_test/classic_B01vsB56.csv",     header = FALSE)
o   <- fread("results/permutation_test/classic_O01vsO20.csv",     header = FALSE)

perm_pvals <- cbind(obo,ob,nbo,nb,o,b)
colnames(perm_pvals) <- c("obo", "ob", "nbo", "nb", "o", "b")

# 3.2 - Adapted ################################################################
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

# Read files and create a new csv with all needed information
obo_adapted <- fread("results/permutation_test/adapted_OBO01vsOBO20.csv", header = FALSE)
ob_adapted  <- fread("results/permutation_test/adapted_OB01vsOB20.csv",   header = FALSE)
nbo_adapted <- fread("results/permutation_test/adapted_nBO01vsnBO56.csv", header = FALSE)
nb_adapted  <- fread("results/permutation_test/adapted_nB01vsnB56.csv",   header = FALSE)
b_adapted   <- fread("results/permutation_test/adapted_B01vsB56.csv",     header = FALSE)
o_adapted   <- fread("results/permutation_test/adapted_O01vsO20.csv",     header = FALSE)

perm_pvals_adapted <-
  cbind(
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

# 4 - Plotting the results =====================================================
cmh_pvals <- readRDS("results/cmh_pvals.rds")
perm_pvals <- fread("results/perm_pvals.csv")

y_limit_up <- 220

width <- 7740
height <- 1440

# 4.1 - Crude ==================================================================
# 4.1.1 - Classic CMH ==========================================================
plot_cmh_classic_OBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20),
    permutation_pvals = perm_pvals$obo,
    title = "Classical CMH test: OBO gen01 vs OBO gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_OBO.png",
  plot = plot_cmh_classic_OBO,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_OB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20),
    permutation_pvals = perm_pvals$ob,
    title = "Classical CMH test: OB gen01 vs OB gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_OB.png",
  plot = plot_cmh_classic_OB,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_nBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_nbo01_vs_nbo56),
    permutation_pvals = perm_pvals$nbo,
    title = "Classical CMH test: nBO gen01 vs nBO gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_nBO.png",
  plot = plot_cmh_classic_nBO,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_nB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_nb01_vs_nb56),
    permutation_pvals = perm_pvals$nb,
    title = "Classical CMH test: nB gen01 vs nB gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_nB.png",
  plot = plot_cmh_classic_nB,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_O <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20),
    permutation_pvals = perm_pvals$o,
    title = "Classical CMH test: O gen01 vs O gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_O.png",
  plot = plot_cmh_classic_O,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_B <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_b01_vs_b56),
    permutation_pvals = perm_pvals$b,
    title = "Classical CMH test: B gen01 vs B gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_B.png",
  plot = plot_cmh_classic_B,
  width = width,
  height = height,
  bg = "white",
  units = "px")

layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_classic_OBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20),
    permutation_pvals = perm_pvals$obo,
    percentage_significance = FALSE,
    title = "Classical CMH test: OBO gen01 vs OBO gen20",
    x_label = FALSE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

grid_plot_cmh_classic_OB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20),
    permutation_pvals = perm_pvals$ob,
    percentage_significance = FALSE,
    title = "Classical CMH test: OB gen01 vs OB gen20",
    x_label = FALSE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

grid_plot_cmh_classic_O <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20),
    permutation_pvals = perm_pvals$o,
    percentage_significance = FALSE,
    title = "Classical CMH test: O gen01 vs O gen20",
    x_label = TRUE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

png(
  filename = "results/figures/cmh_crude/classic/cmh_classic_OBO_OB_O_piled.png",
  width = 1800,
  height = 900)

grid.arrange(
  grid_plot_cmh_classic_OBO,
  grid_plot_cmh_classic_OB,
  grid_plot_cmh_classic_O,
  layout_matrix = layout)

dev.off()

# 4.1.2 - Adapted CMH ==========================================================
plot_cmh_adapted_OBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: OBO gen01 vs OBO gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_OBO.png",
  plot = plot_cmh_adapted_OBO,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_OB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: OB gen01 vs OB gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_OB.png",
  plot = plot_cmh_adapted_OB,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_nBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_nbo01_vs_nbo56),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: nBO gen01 vs nBO gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_nBO.png",
  plot = plot_cmh_adapted_nBO,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_nB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_nb01_vs_nb56),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: nB gen01 vs nB gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_nB.png",
  plot = plot_cmh_adapted_nB,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_O <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: O gen01 vs O gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_O.png",
  plot = plot_cmh_adapted_O,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_B <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_b01_vs_b56),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: B gen01 vs B gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_B.png",
  plot = plot_cmh_adapted_B,
  width = width,
  height = height,
  bg = "white",
  units = "px")

layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_adapted_OBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: OBO gen01 vs OBO gen20",
    x_label = FALSE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

grid_plot_cmh_adapted_OB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: OB gen01 vs OB gen20",
    x_label = FALSE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

grid_plot_cmh_adapted_O <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: O gen01 vs O gen20",
    x_label = TRUE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

png(
  filename = "results/figures/cmh_crude/adapted/cmh_adapted_OBO_OB_O_piled.png",
  width = 1800,
  height = 900)

grid.arrange(
  grid_plot_cmh_adapted_OBO,
  grid_plot_cmh_adapted_OB,
  grid_plot_cmh_adapted_O,
  layout_matrix = layout)

dev.off()

# 4.2 - Scaled data ============================================================
# 4.2.1 - Classic CMH Scaled ===================================================
plot_cmh_classic_OBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: OBO gen01 vs OBO gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_OBO_scaled.png",
       plot = plot_cmh_classic_OBO_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_OB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: OB gen01 vs OB gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_OB_scaled.png",
       plot = plot_cmh_classic_OB_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_nBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_nbo01_vs_nbo56_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: nBO gen01 vs nBO gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_nBO_scaled.png",
       plot = plot_cmh_classic_nBO_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_nB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_nb01_vs_nb56_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: nB gen01 vs nB gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_nB_scaled.png",
       plot = plot_cmh_classic_nB_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_O_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_O_scaled.png",
       plot = plot_cmh_classic_O_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_B_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_b01_vs_b56_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: B gen01 vs B gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_B_scaled.png",
       plot = plot_cmh_classic_B_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

# Grid arrange  
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_classic_OBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: OBO gen01 vs OBO gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_classic_OB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: OB gen01 vs OB gen20",
           x_label = FALSE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_classic_O_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_scaled/classic/cmh_classic_OBO_OB_O_scaled_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_classic_OBO_scaled,
             grid_plot_cmh_classic_OB_scaled,
             grid_plot_cmh_classic_O_scaled,
             layout_matrix = layout)

dev.off()

# 4.2.2 - Adapted CMH Scaled ===================================================
plot_cmh_adapted_OBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: OBO gen01 vs OBO gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_OBO_scaled.png",
       plot = plot_cmh_adapted_OBO_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_OB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: OB gen01 vs OB gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_OB_scaled.png",
       plot = plot_cmh_adapted_OB_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_nBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_nbo01_vs_nbo56_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: nBO gen01 vs nBO gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_nBO_scaled.png",
       plot = plot_cmh_adapted_nBO_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_nB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_nb01_vs_nb56_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: nB gen01 vs nB gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_nB_scaled.png",
       plot = plot_cmh_adapted_nB_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_O_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_O_scaled.png",
       plot = plot_cmh_adapted_O_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_B_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_b01_vs_b56_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: B gen01 vs B gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_B_scaled.png",
       plot = plot_cmh_adapted_B_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_adapted_OBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: OBO gen01 vs OBO gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_adapted_OB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: OB gen01 vs OB gen20",
           x_label = FALSE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_adapted_O_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_scaled/adapted/cmh_adapted_OBO_OB_O_scaled_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_adapted_OBO_scaled,
             grid_plot_cmh_adapted_OB_scaled,
             grid_plot_cmh_adapted_O_scaled,
             layout_matrix = layout)

dev.off()

# 4.3 - FDR Corrected ==========================================================
# 4.3.1 - Classic CMH FDR ======================================================
plot_cmh_classic_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_obo01_vs_obo20, method = "BH")),
                   permutation_pvals = perm_pvals$obo,
                   title = "Classical CMH test, FDR corrected: OBO gen01 vs OBO gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_OBO.png",
       plot = plot_cmh_classic_OBO,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_ob01_vs_ob20, method = "BH")),
                   permutation_pvals = perm_pvals$ob,
                   title = "Classical CMH test, FDR corrected: OB gen01 vs OB gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_OB.png",
       plot = plot_cmh_classic_OB,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_nBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_nbo01_vs_nbo56, method = "BH")),
                   permutation_pvals = perm_pvals$nbo,
                   title = "Classical CMH test, FDR corrected: nBO gen01 vs nBO gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_nBO.png",
       plot = plot_cmh_classic_nBO,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_nB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_nb01_vs_nb56, method = "BH")),
                   permutation_pvals = perm_pvals$nb,
                   title = "Classical CMH test, FDR corrected: nB gen01 vs nB gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_nB.png",
       plot = plot_cmh_classic_nB,
       width = width,
       height = height,
       bg = "white",
       units = "px")


plot_cmh_classic_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_o01_vs_o20, method = "BH")),
                   permutation_pvals = perm_pvals$o,
                   title = "Classical CMH test, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_O.png",
       plot = plot_cmh_classic_O,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_B <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_b01_vs_b56, method = "BH")),
                   permutation_pvals = perm_pvals$b,
                   title = "Classical CMH test, FDR corrected: B gen01 vs B gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_B.png",
       plot = plot_cmh_classic_B,
       width = width,
       height = height,
       bg = "white",
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_classic_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_obo01_vs_obo20, method = "BH")),
                   permutation_pvals = perm_pvals$obo,
                   title = "Classical CMH test, FDR corrected: OBO gen01 vs OBO gen20",
                   x_label = FALSE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_classic_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_ob01_vs_ob20, method = "BH")),
                   permutation_pvals = perm_pvals$ob,
                   title = "Classical CMH test, FDR corrected: OB gen01 vs OB gen20",
                   x_label = FALSE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_classic_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_o01_vs_o20, method = "BH")),
                   permutation_pvals = NULL,
                   title = "Classical CMH test, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_fdr/classic/cmh_classic_OBO_OB_O_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_classic_OBO,
             grid_plot_cmh_classic_OB,
             grid_plot_cmh_classic_O,
             layout_matrix = layout)

dev.off()

# 4.3.2 - Adapted CMH FDR ======================================================
plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_obo01_vs_obo20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: OBO gen01 vs OBO gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_OBO.png",
       plot = plot_cmh_adapted_OBO,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_ob01_vs_ob20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: OB gen01 vs OB gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_OB.png",
       plot = plot_cmh_adapted_OB,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_nBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_nbo01_vs_nbo56, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: nBO gen01 vs nBO gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_nBO.png",
       plot = plot_cmh_adapted_nBO,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_nB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_nb01_vs_nb56, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: nB gen01 vs nB gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_nB.png",
       plot = plot_cmh_adapted_nB,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_o01_vs_o20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_O.png",
       plot = plot_cmh_adapted_O,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_B <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_b01_vs_b56, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: B gen01 vs B gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_B.png",
       plot = plot_cmh_adapted_B,
       width = width,
       height = height,
       bg = "white",
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_obo01_vs_obo20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: OBO gen01 vs OBO gen20",
                   x_label = FALSE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_ob01_vs_ob20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: OB gen01 vs OB gen20",
                   x_label = FALSE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_o01_vs_o20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_fdr/adapted/cmh_adapted_OBO_OB_O_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_adapted_OBO,
             grid_plot_cmh_adapted_OB,
             grid_plot_cmh_adapted_O,
             layout_matrix = layout)

dev.off()

# 4.4 - Scaled and FDR =========================================================
# Pending all plots!!!
# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)
y_limit_up <- 220

grid_plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_obo01_vs_obo20_scaled, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, scaled, FDR corrected: OBO gen01 vs OBO gen20",
                   x_label = FALSE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_ob01_vs_ob20_scaled, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, scaled, FDR corrected: OB gen01 vs OB gen20",
                   x_label = FALSE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_o01_vs_o20_scaled, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, scaled, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_fdr_scaled/adapted/cmh_adapted_OBO_OB_O_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_adapted_OBO,
             grid_plot_cmh_adapted_OB,
             grid_plot_cmh_adapted_O,
             layout_matrix = layout)

dev.off()


# 4.4 - QQ Plots ================================================================

n <- nrow(cmh_pvals)

# Create a vector of N values evenly spaces from 1 to 1 / N
N <- sort(-log10(seq(1, n) / n))


p <- sort(-log10(cmh_pvals$cmh_adapted_o01_vs_o20_scaled))

# 5 - PCA Analysis =============================================================
snp_table_shahrestani <- 
  as.data.frame(fread("data/processed/processed_snps_abcd_shahrestani.csv"))

pca_data <- PreparePca(snp_table_shahrestani)

pca_plot_labeled <- PlotPca(pca_data, label = TRUE)
ggsave("results/figures/pca_plot_labeled.png",
       plot = pca_plot_labeled,
       width = 2000,
       height = 2000,
       units = "px")

pca_plot_unlabeled <- PlotPca(pca_data, label = FALSE)
ggsave("results/figures/pca_plot_unlabeled.png",
       plot = pca_plot_unlabeled,
       width = width,
       height = 1600,
       units = "px")

# 6 - Frequency Analysis =======================================================
# Load files
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")

# Check freqs for all SNPs
freq <- GetFreq(snp_table_shahrestani)

# Create a layout with two rows and 5 columns
layout <- matrix(c(1,2,3,4,5,6,7,8,9,10), ncol = 5, byrow = TRUE)

# Plot OBO freqs
obo_rep01_gen01_freq <- 
  ggplot(freq, aes(alt_OBO_rep01_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep02_gen01_freq <-
  ggplot(freq, aes(alt_OBO_rep02_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep03_gen01_freq <-
  ggplot(freq, aes(alt_OBO_rep03_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep04_gen01_freq <-
  ggplot(freq, aes(alt_OBO_rep04_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep05_gen01_freq <-
  ggplot(freq, aes(alt_OBO_rep05_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep01_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep01_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep02_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep02_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep03_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep03_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep04_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep04_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep05_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep05_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

png(filename = "results/figures/freq_hist/obo_freq_hist_all.png",
    width = 1800,
    height = 900)

grid.arrange(obo_rep01_gen01_freq,
             obo_rep02_gen01_freq,
             obo_rep03_gen01_freq,
             obo_rep04_gen01_freq,
             obo_rep05_gen01_freq,
             obo_rep01_gen20_freq,
             obo_rep02_gen20_freq,
             obo_rep03_gen20_freq,
             obo_rep04_gen20_freq,
             obo_rep05_gen20_freq,
             layout_matrix = layout)

dev.off()

# Plot OB freqs
OB_rep01_gen01_freq <- 
  ggplot(freq, aes(alt_OB_rep01_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep02_gen01_freq <-
  ggplot(freq, aes(alt_OB_rep02_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep03_gen01_freq <-
  ggplot(freq, aes(alt_OB_rep03_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep04_gen01_freq <-
  ggplot(freq, aes(alt_OB_rep04_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep05_gen01_freq <-
  ggplot(freq, aes(alt_OB_rep05_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep01_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep01_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep02_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep02_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep03_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep03_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep04_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep04_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep05_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep05_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

png(filename = "results/figures/freq_hist/OB_freq_hist_all.png",
    width = 1800,
    height = 900)

grid.arrange(OB_rep01_gen01_freq,
             OB_rep02_gen01_freq,
             OB_rep03_gen01_freq,
             OB_rep04_gen01_freq,
             OB_rep05_gen01_freq,
             OB_rep01_gen20_freq,
             OB_rep02_gen20_freq,
             OB_rep03_gen20_freq,
             OB_rep04_gen20_freq,
             OB_rep05_gen20_freq,
             layout_matrix = layout)

dev.off()

# Now do the same thing for only significant SNPs
# Check freqs for all SNPs
# Add the CMH vals to the snp_table so we can filter them all together
cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_shahrestani <- cbind(snp_table_shahrestani, cmh_pvals)

threshold <- 1e-100

filtered_snp_table <- snp_table_shahrestani[which(snp_table_shahrestani$cmh_adapted_o01_vs_o20 < threshold),]

freq_significant <- GetFreq(filtered_snp_table)

# Create a layout with two rows and 5 columns
layout <- matrix(c(1,2,3,4,5,6,7,8,9,10), ncol = 5, byrow = TRUE)

# Plot OBO freqs
obo_rep01_gen01_freq_significant <- 
  ggplot(freq_significant, aes(alt_OBO_rep01_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep02_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep02_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep03_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep03_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep04_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep04_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep05_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep05_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep01_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep01_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep02_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep02_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep03_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep03_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep04_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep04_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep05_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep05_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

png(filename = "results/figures/freq_hist/obo_freq_significant_hist.png",
    width = 1800,
    height = 900)

grid.arrange(obo_rep01_gen01_freq_significant,
             obo_rep02_gen01_freq_significant,
             obo_rep03_gen01_freq_significant,
             obo_rep04_gen01_freq_significant,
             obo_rep05_gen01_freq_significant,
             obo_rep01_gen20_freq_significant,
             obo_rep02_gen20_freq_significant,
             obo_rep03_gen20_freq_significant,
             obo_rep04_gen20_freq_significant,
             obo_rep05_gen20_freq_significant,
             layout_matrix = layout)

dev.off()

# Plot OB freqs
OB_rep01_gen01_freq_significant <- 
  ggplot(freq_significant, aes(alt_OB_rep01_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep02_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep02_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep03_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep03_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep04_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep04_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep05_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep05_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep01_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep01_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep02_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep02_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep03_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep03_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep04_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep04_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep05_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep05_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

png(filename = "results/figures/freq_hist/OB_freq_significant_hist.png",
    width = 1800,
    height = 900)

grid.arrange(OB_rep01_gen01_freq_significant,
             OB_rep02_gen01_freq_significant,
             OB_rep03_gen01_freq_significant,
             OB_rep04_gen01_freq_significant,
             OB_rep05_gen01_freq_significant,
             OB_rep01_gen20_freq_significant,
             OB_rep02_gen20_freq_significant,
             OB_rep03_gen20_freq_significant,
             OB_rep04_gen20_freq_significant,
             OB_rep05_gen20_freq_significant,
             layout_matrix = layout)

dev.off()


# 7 - Allele Trajectory Analysis ===============================================

# Data loading
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")


PhaseSnps <- function(freq) {
  
  print("Warning: phasing SNPs put them on a different col order")
  print("Now, we have all gen01 columns and then all gen20 columns")
  
  phased_snp_table_gen01 <- freq[,grep("gen01", colnames(freq))]
  phased_snp_table_gen20 <- freq[,grep("gen20|gen56", colnames(freq))]
  
  marked_for_phasing <- phased_snp_table_gen01 > 0.5
  
  # If marked for phasing, subtract 1 from the value
  phased_snp_table_gen01[marked_for_phasing] <- 
    1 - phased_snp_table_gen01[marked_for_phasing]
  
  phased_snp_table_gen20[marked_for_phasing] <- 
    1 - phased_snp_table_gen20[marked_for_phasing]
  
  phased_snp_table <- cbind(phased_snp_table_gen01, phased_snp_table_gen20)
  
  return(phased_snp_table)
}

cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_shahrestani <- cbind(snp_table_shahrestani, cmh_pvals)

threshold <- 1e-100

filtered_snp_table <- snp_table_shahrestani[which(snp_table_shahrestani$cmh_adapted_o01_vs_o20 < threshold),]
filtered_snp_table <- dplyr::select(filtered_snp_table, -contains("nB"))

freq_significant <- GetFreq(filtered_snp_table)

phased_freq_significant <- PhaseSnps(freq_significant)

#phased_freq_significant_OBO <- phased_freq_significant[,grep("OBO", colnames(phased_freq_significant))]
#phased_freq_significant_OB <- phased_freq_significant[,grep("OB_", colnames(phased_freq_significant))]

plots <- list()
for (i in 0:9) {
  # Define y axis
  rep_i_01 <- phased_freq_significant[1+i]
  print(colnames(rep_i_01))
  rep_i_20 <- phased_freq_significant[11+i]
  print(colnames(rep_i_20))
  
  colnames(rep_i_01) <- "freq"
  colnames(rep_i_20) <- "freq"
  
  plotting_df <- cbind(rep_i_01, rep_i_20)
  
  # png(filename = paste("results/figures/allele_trajectory/rep", i+1, ".png", sep = ""),
  #     width = 1800,
  #     height = 900)
  
  plot(1,
       type = "n",
       xlab = "", 
       ylab = "",
       xlim = c(1, 20),  
       ylim = c(0, 1),
       xaxt = "n",
       main = paste("rep", i+1),
       cex.axis = 2
  )
  for (row in 1:nrow(plotting_df)) {
    lines(x = c(1,20), 
          y = plotting_df[row,])
  }
  
  p <- recordPlot()
  
  # dev.off()
  
  plots[[i+1]] <- p
  print(i+1)
}


# 7.1 - Delta Statistic ========================================================
# Data loading
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

freq <- GetFreq(snp_table_shahrestani)
freq <- dplyr::select(freq, -contains("nB"))

phased_freq <- PhaseSnps(freq)

# Create a delta dataframe, which will calculate the delta in frequency for all the replicates
delta_df <- data.frame(matrix(NA, nrow = nrow(freq), ncol = 10))

gen01 <- phased_freq[,grep("gen01", colnames(phased_freq))]
gen20 <- phased_freq[,grep("gen20|gen56", colnames(phased_freq))]

for (i in 1:10) {
  delta_df[,i] <- gen20[,i] - gen01[,i]
}

colnames(delta_df) <- colnames(freq[1:10])
delta_df <- abs(delta_df)

cmh_pvals <- readRDS("results/cmh_pvals.rds")
cmh_pvals$delta <- rowSums(delta_df)/ncol(delta_df)

y_limit_up <- 220
grid_plot_cmh_adapted_O_scaled_fdr <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(p.adjust(cmh_pvals$cmh_adapted_o01_vs_o20_scaled, method = "BH")),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test, scaled, FDR corrected: O gen01 vs O gen20, colored for delta > 0.40",
    x_label = TRUE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)


# Lets color the top 10% higher delta values in red
grid_plot_cmh_adapted_O_scaled_fdr <-
  grid_plot_cmh_adapted_O_scaled_fdr +
  aes(color = delta > 0.40) +  # Add color mapping
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(color = "Higher Delta")  # Update legend label


grid_plot_cmh_adapted_O_scaled_fdr

# ==============================================================================
  
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")

# Add the CMH vals to the snp_table so we can filter them all together
cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_shahrestani <- cbind(snp_table_shahrestani, cmh_pvals)

threshold <- 1e-100

filtered_snp_table <- snp_table_shahrestani[which(snp_table_shahrestani$cmh_adapted_o01_vs_o20 < threshold),]

freq <- GetFreq(filtered_snp_table)

# Problem: most significant SNPs start from a fixed frequency all O-type populations
# If that were true, it would mean that the exact same mutation happened in all experimental populations, which is extremely unlikely
# This is probably an artifact of our data processing workflow

# out of the significant snps, make a histogram of the 

# This will show you that the most significant SNPs start from a fixed frequency
snp_table_significant <- 
  snp_table_shahrestani[which(snp_table_shahrestani$cmh_adapted_o01_vs_o20 < threshold),]
DiagnoseSnps(snp_table_significant)

# Let's filter out the SNPs that start from a fixed frequency
filtered_snp_table <- FilterOutFixedSnps(snp_table_shahrestani)

# Let's see how a manhattan plot with only these SNPs looks like
# The perm_pval can't be used here, as the permutation test was run with all SNPs
# I would need to run a new permutation test with only these 407,678 SNPs to be sure.
filtered_manhplot <-
  GetManhattanPlot(
    my_dataframe = filtered_snp_table,
    Y = -log10(filtered_snp_table$cmh_adapted_o01_vs_o20),
    #permutation_pvals = perm_pvals$o,
    title = "Filtered Manhattan plot - O gen01 vs O gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = 100,
    y_limit_down = 0
  )

# Lets paint the SNPs that are above the threshold in red
filtered_manhplot <- filtered_manhplot +
  aes(color = -log10(cmh_adapted_o01_vs_o20) > -log10(threshold)) +  # Add color mapping
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(color = "Above Threshold")  # Update legend label

filtered_manhplot

# Let's see the how the most significant SNPs look like now
filtered_snp_table_significant <-
  filtered_snp_table[which(filtered_snp_table$cmh_adapted_o01_vs_o20 < threshold),]

DiagnoseSnps(filtered_snp_table_significant)

# Now we have interesting SNPs to look at. Let's do an enrichment analysis
GO_dataframe <- filtered_snp_table_significant[c("CHROM", "POS", "cmh_adapted_o01_vs_o20")]
GO_dataframe$cmh_adapted_o01_vs_o20 <- -log10(GO_dataframe$cmh_adapted_o01_vs_o20)
GO_dataframe$coordinate <- paste(GO_dataframe$CHROM, ":", GO_dataframe$POS, sep = "")

# I've manually identified peaks
peaks <- c(
  "peak_1 " = "2L:2901919",
  "peak_2 " = "2L:9348851",
  "peak_3 " = "2L:10187059",
  "peak_4 " = "2L:17541874",
  "peak_5 " = "2L:21931301",
  "peak_6 " = "2R:4559132",
  "peak_7 " = "2R:8553887",
  "peak_8 " = "2R:9814502",
  "peak_9 " = "2R:12964273",
  "peak_10" = "2R:17126389", # Check this one
  "peak_11" = "2R:21846663", # This is a weird dot that falls on the line
  "peak_12" = "2R:25001106",  # This one is also weird
  "peak_13" = "3L:9338136",
  "peak_14" = "3L:9787463",
  "peak_15" = "3L:11211200",
  "peak_16" = "3L:13404028",
  "peak_17" = "3L:15767228",
  "peak_18" = "3L:16327479",
  "peak_19" = "3L:19958564",
  "peak_20" = "3L:22980227",
  "peak_21" = "3L:24755959",
  "peak_22" = "3R:4343926",
  "peak_23" = "3R:5218031",
  "peak_24" = "3R:5231920",
  "peak_25" = "3R:5860717",
  "peak_26" = "3R:6977845",
  "peak_27" = "3R:7220441",
  "peak_28" = "3R:10865062",
  "peak_29" = "3R:12723625",
  "peak_30" = "3R:13064936",
  "peak_31" = "3R:15178836",
  "peak_32" = "3R:17175981",
  "peak_33" = "3R:18127342",
  "peak_34" = "3R:23236137",
  "peak_35" = "3R:27000514",
  "peak_36" = "3R:28166423",
  "peak_37" = "3R:28989889",
  "peak_38" = "3R:30694517",
  "peak_39" = "X:2499092",
  "peak_40" = "X:3532997",
  "peak_41" = "X:5518016",
  "peak_42" = "X:6194353"
)

# Now we can use biomaRt to get a gene list for those regions.
# We can later use that genelist in a website like GOrilla and see if there is any enrichment
ensembl <- useEnsembl(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl")
