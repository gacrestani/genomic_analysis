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

# 1 - Process raw data =========================================================
snp_table_shahrestani <- as.data.frame(ReadAndPrepare(mode = "shahrestani"))
saveRDS(snp_table_shahrestani,
        "data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <- as.data.frame(ReadAndPrepare(mode = "regimes"))
saveRDS(snp_table_regimes,
        "data/processed/processed_snps_abcd_regimes.rds")

# 2 - CMH tests ================================================================
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

# Initiates an empty dataframe with the same number of rows as our snp tables
cmh_pvals <- data.frame(matrix(NA,
                               nrow = nrow(snp_table_shahrestani),
                               ncol = 0))

cmh_pvals$ABS_POS <- snp_table_shahrestani$ABS_POS
cmh_pvals$CHROM <- snp_table_shahrestani$CHROM

# Classical CMH test
# NOTE: These functions take a while to run (even running in parallel)
cmh_pvals$cmh_classic_obo01_vs_obo20 <-
  ClassicalCmhTest(snp_table = snp_table_shahrestani,
                   treatment1 = "OBO",
                   gen1 = "01",
                   treatment2 = "OBO",
                   gen2 = "20")

cmh_pvals$cmh_classic_ob01_vs_ob20 <-
  ClassicalCmhTest(snp_table = snp_table_shahrestani,
                   treatment1 = "OB",
                   gen1 = "01",
                   treatment2 = "OB",
                   gen2 = "20")

cmh_pvals$cmh_classic_nbo01_vs_nbo56 <-
  ClassicalCmhTest(snp_table = snp_table_shahrestani,
                   treatment1 = "nBO",
                   gen1 = "01",
                   treatment2 = "nBO",
                   gen2 = "56")

cmh_pvals$cmh_classic_nb01_vs_nb56 <-
  ClassicalCmhTest(snp_table = snp_table_shahrestani,
                   treatment1 = "nB",
                   gen1 = "01",
                   treatment2 = "nB",
                   gen2 = "56")

cmh_pvals$cmh_classic_o01_vs_o20 <-
  ClassicalCmhTest(snp_table = snp_table_regimes,
                   treatment1 = "O",
                   gen1 = "01",
                   treatment2 = "O",
                   gen2 = "20")

cmh_pvals$cmh_classic_b01_vs_b56 <-
  ClassicalCmhTest(snp_table = snp_table_regimes,
                   treatment1 = "B",
                   gen1 = "01",
                   treatment2 = "B",
                   gen2 = "56")

# Adapted CMH test
# These are faster than the classical tests
cmh_pvals$cmh_adapted_obo01_vs_obo20 <-
  AdaptedCmhTest(snp_table = snp_table_shahrestani,
                 treatment1 = "OBO",
                 gen1 = "01",
                 treatment2 = "OBO",
                 gen2 = "20")

cmh_pvals$cmh_adapted_ob01_vs_ob20 <-
  AdaptedCmhTest(snp_table = snp_table_shahrestani,
                 treatment1 = "OB",
                 gen1 = "01",
                 treatment2 = "OB",
                 gen2 = "20")

cmh_pvals$cmh_adapted_nbo01_vs_nbo56 <-
  AdaptedCmhTest(snp_table = snp_table_shahrestani,
                 treatment1 = "nBO",
                 gen1 = "01",
                 treatment2 = "nBO",
                 gen2 = "56")

cmh_pvals$cmh_adapted_nb01_vs_nb56 <-
  AdaptedCmhTest(snp_table = snp_table_shahrestani,
                 treatment1 = "nB",
                 gen1 = "01",
                 treatment2 = "nB",
                 gen2 = "56")

cmh_pvals$cmh_adapted_o01_vs_o20 <-
  AdaptedCmhTest(snp_table = snp_table_regimes,
                 treatment1 = "O",
                 gen1 = "01",
                 treatment2 = "O",
                 gen2 = "20")

cmh_pvals$cmh_adapted_b01_vs_b56 <-
  AdaptedCmhTest(snp_table = snp_table_regimes,
                 treatment1 = "B",
                 gen1 = "01",
                 treatment2 = "B",
                 gen2 = "56")

saveRDS(cmh_pvals, "results/cmh_pvals.rds")

# I have p-values so small that the multiple-test correction methods used by
# other papers are not suitable. If I apply the q-value method and consider
# everything < 0.05 as significant, the whole genome will have significant snps.
# Possible problems: my Ne estimationsa are VERY low, like 30. This does not
# represent reality. I should check in with Molly about this

# # Q-value experiment
# q_values <- qvalue(cmh_pvals$cmh_classic_o01_vs_o20)
# 
# plot_cmh_classic_O_q <-
#   GetManhattanPlot(my_dataframe = cmh_pvals,
#                    Y = -log10(q_values$qvalue),
#                    permutation_pvals = perm_pvals$o,
#                    title = "Classical CMH test: O gen01 vs O gen20 - q-value",
#                    x_label = TRUE,
#                    y_label = NULL,
#                    palette = "blue",
#                    y_limit_up = y_limit_up,
#                    y_limit_down = 0)
# ggsave("results/figures/cmh_classic_O_q.png",
#        plot = plot_cmh_classic_O_q,
#        width = 1600,
#        height = 1200,
#        bg = "white",
#        units = "px")
# 
# sig <- filter(q_values$qvalues, qval < 0.05)
# 
# q_vals <- q_values$qvalues
# sig <- q_vals[q_vals < 0.05]
# 
# q_values_adapted <- qvalue(cmh_pvals$cmh_adapted_o01_vs_o20)
# plot(q_values_adapted)
# 
# plot(x = cmh_pvals$cmh_adapted_o01_vs_o20, y = q_values_adapted$qvalues)
# 
# plot(x = cmh_pvals$ABS_POS, y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20))
# abline(h = -log10(0.05/length(cmh_pvals$cmh_classic_o01_vs_o20)), col = "red")
# 
# plot(x = cmh_pvals$ABS_POS, y = -log10(q_values_adapted$qvalue))
# abline(h = -log10(0.05), col = "red")
# 
# q001 <- quantile(as.numeric(cmh_pvals$cmh_adapted_o01_vs_o20), probs = 0.001)
# plot(x = cmh_pvals$ABS_POS, y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20))
# abline(h = -log10(q001), col = "red")

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

fwrite(perm_pvals, file = "results/perm_pvals.csv")

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
  gen2 = "20"
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
  gen2 = "20"
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
  gen2 = "56"
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
  gen2 = "56"
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
  gen2 = "20"
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
  gen2 = "56"
)

# Read files and create a new csv with all needed information
obo <- fread("results/permutation_test/adapted_OBO01vsOBO20.csv", header = FALSE)
ob  <- fread("results/permutation_test/adapted_OB01vsOB20.csv",   header = FALSE)
nbo <- fread("results/permutation_test/adapted_nBO01vsnBO56.csv", header = FALSE)
nb  <- fread("results/permutation_test/adapted_nB01vsnB56.csv",   header = FALSE)
b   <- fread("results/permutation_test/adapted_B01vsB56.csv",     header = FALSE)
o   <- fread("results/permutation_test/adapted_O01vsO20.csv",     header = FALSE)

perm_pvals <- cbind(obo,ob,nbo,nb,o,b)
colnames(perm_pvals) <- c("obo", "ob", "nbo", "nb", "o", "b")

fwrite(perm_pvals, file = "results/perm_pvals.csv")




# 4 - Plotting the results =====================================================
cmh_pvals <- readRDS("results/cmh_pvals.rds")
perm_pvals <- fread("results/perm_pvals.csv")

y_limit_up <- 200

# 4.1 - Classic CMH ====
plot_cmh_classic_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20),
           permutation_pvals = perm_pvals$obo,
           title = "Classical CMH test: OBO gen01 vs OBO gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_classic_OBO.png",
       plot = plot_cmh_classic_OBO,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_classic_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20),
           permutation_pvals = perm_pvals$ob,
           title = "Classical CMH test: OB gen01 vs OB gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_classic_OB.png",
       plot = plot_cmh_classic_OB,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_classic_nBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_nbo01_vs_nbo56),
           permutation_pvals = perm_pvals$nbo,
           title = "Classical CMH test: nBO gen01 vs nBO gen56",
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_classic_nBO.png",
       plot = plot_cmh_classic_nBO,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_classic_nB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_nb01_vs_nb56),
           permutation_pvals = perm_pvals$nb,
           title = "Classical CMH test: nB gen01 vs nB gen56",
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_classic_nB.png",
       plot = plot_cmh_classic_nB,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")


plot_cmh_classic_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20),
           permutation_pvals = perm_pvals$o,
           title = "Classical CMH test: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_classic_O.png",
       plot = plot_cmh_classic_O,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_classic_B <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_b01_vs_b56),
           permutation_pvals = perm_pvals$b,
           title = "Classical CMH test: B gen01 vs B gen56",
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_classic_B.png",
       plot = plot_cmh_classic_B,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_classic_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20),
           permutation_pvals = perm_pvals$obo,
           title = "Classical CMH test: OBO gen01 vs OBO gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_classic_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20),
           permutation_pvals = perm_pvals$ob,
           title = "Classical CMH test: OB gen01 vs OB gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_classic_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20),
           permutation_pvals = perm_pvals$o,
           title = "Classical CMH test: O gen01 vs O gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_classic_OBO_OB_O_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_classic_OBO,
             grid_plot_cmh_classic_OB,
             grid_plot_cmh_classic_O,
             layout_matrix = layout)

dev.off()

# 4.2 - Adapted CMH ====
plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: OBO gen01 vs OBO gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_OBO.png",
       plot = plot_cmh_adapted_OBO,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: OB gen01 vs OB gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_OB.png",
       plot = plot_cmh_adapted_OB,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_adapted_nBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_nbo01_vs_nbo56),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: nBO gen01 vs nBO gen56",
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_nBO.png",
       plot = plot_cmh_adapted_nBO,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_adapted_nB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_nb01_vs_nb56),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: nB gen01 vs nB gen56",
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_nB.png",
       plot = plot_cmh_adapted_nB,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_O.png",
       plot = plot_cmh_adapted_O,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

plot_cmh_adapted_B <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_b01_vs_b56),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: B gen01 vs B gen56",
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_B.png",
       plot = plot_cmh_adapted_B,
       width = 1600,
       height = 1200,
       bg = "white",
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: OBO gen01 vs OBO gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: OB gen01 vs OB gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
           permutation_pvals = NULL,
           percentage_significance = 0.001,
           title = "Adapted CMH test: O gen01 vs O gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_adapted_OBO_OB_O_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_adapted_OBO,
             grid_plot_cmh_adapted_OB,
             grid_plot_cmh_adapted_O,
             layout_matrix = layout)

dev.off()

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
       width = 1600,
       height = 1600,
       units = "px")

# 6 - Gene Enrichment Analysis =================================================
# Load files
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")

# Add the CMH vals to the snp_table so we can filter them all together
cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_shahrestani <- cbind(snp_table_shahrestani, cmh_pvals)

# Let's create a threshold of significance
threshold <- quantile(
  as.numeric(filtered_snp_table$cmh_adapted_o01_vs_o20),
  probs = 0.001
)

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
