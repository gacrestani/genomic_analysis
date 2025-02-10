# ==============================================================================
# This script is the main script that runs the entire analysis. It reads the raw
# snp_table, processes it, and runs the CMH tests. The results are saved in the
# results folder.
# ==============================================================================

# Read raw snp_table and create a processed version ============================
source("scripts/01_data_preparation.R")

snp_table_shahrestani <- as.data.frame(ReadAndPrepare(mode = "shahrestani"))
fwrite(snp_table_shahrestani,
       file = "data/processed/processed_snps_abcd_shahrestani.csv",
       sep = ",",
       row.names = FALSE)

snp_table_regimes <- as.data.frame(ReadAndPrepare(mode = "regimes"))
fwrite(snp_table_regimes,
       file = "data/processed/processed_snps_abcd_regimes.csv",
       sep = ",",
       row.names = FALSE)
# ==============================================================================

# Running the CMH tests ========================================================
source("scripts/02_cmh_tests.R")

snp_table_shahrestani <- 
  as.data.frame(fread("data/processed/processed_snps_abcd_shahrestani.csv"))

snp_table_regimes <-
  as.data.frame(fread("data/processed/processed_snps_abcd_regimes.csv"))

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

fwrite(cmh_pvals,
       "results/cmh_pvals.csv",
       sep = ",",
       row.names = FALSE)

# Plotting the results =========================================================
source("scripts/03_plot_functions.R")

cmh_pvals <- fread("results/cmh_pvals.csv")
perm_pvals <- fread("results/perm_pvals.csv")

y_limit_up <- 200

# Classic CMH ====
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

# Adapted CMH ====
plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20),
           permutation_pvals = NULL,
           title = NULL,
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_OBO.png",
       plot = plot_cmh_adapted_OBO,
       width = 1600,
       height = 1200,
       units = "px")

plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20),
           permutation_pvals = NULL,
           title = NULL,
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_OB.png",
       plot = plot_cmh_adapted_OB,
       width = 1600,
       height = 1200,
       units = "px")

plot_cmh_adapted_nBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_nbo01_vs_nbo56),
           permutation_pvals = NULL,
           title = NULL,
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_nBO.png",
       plot = plot_cmh_adapted_nBO,
       width = 1600,
       height = 1200,
       units = "px")

plot_cmh_adapted_nB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_nb01_vs_nb56),
           permutation_pvals = NULL,
           title = NULL,
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_nB.png",
       plot = plot_cmh_adapted_nB,
       width = 1600,
       height = 1200,
       units = "px")

plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
           permutation_pvals = NULL,
           title = NULL,
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_O.png",
       plot = plot_cmh_adapted_O,
       width = 1600,
       height = 1200,
       units = "px")

plot_cmh_adapted_B <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_b01_vs_b56),
           permutation_pvals = NULL,
           title = NULL,
           x_label = TRUE,
           y_label = NULL,
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_adapted_B.png",
       plot = plot_cmh_adapted_B,
       width = 1600,
       height = 1200,
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20),
           permutation_pvals = NULL,
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

# PCA Analysis =================================================================
source("scripts/02_cmh_tests.R")
source("scripts/04_pca_analysis.R")

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
