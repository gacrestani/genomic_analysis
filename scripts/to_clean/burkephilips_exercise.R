# Lib imports ####
library(ACER)
library(poolSeq)
library(tidyverse)
library(gridExtra)

# Functions
getAdaptedCMH <- function(data, replicates) {
  
  # Set Replicates for debug run
  #replicates <- control_replicates
  
  # Creates two datasets - allele frequency and coverage
  freq <- data[grep("alt_", colnames(data))]
  cov <- data[grep("N_", colnames(data))]
  
  # Gets Ne (P.planII) in a list, jumping by steps of 3 and therefore accounting for all samples
  Ne <- c()
  for (i in seq(1, ncol(cov), by=3)) {
    unlisted_ne <- estimateNe(p0 = freq[,i],
                              pt = freq[,(i+2)],
                              cov0 = cov[,i],
                              covt = cov[,(i+2)],
                              t = 200,
                              method=c("P.planII"),
                              Ncensus=NA,
                              poolSize=c(1e+7, 1e+7))
    
    Ne <- c(Ne, as.integer(unname(unlisted_ne)))
  }
  
  # Gets the p_vals vector using the just-calculated Ne and accounting for the chosen replicates
  # This will not work on the manual run because data is not filtered correctly.
  # However, this does run when being called by the getCmhPlot, as there is a filtering step there
  p_vals <- adapted.cmh.test(freq = as.matrix(freq),
                             coverage = as.matrix(cov),
                             Ne = Ne,
                             gen = c(1,100,200),
                             repl = as.numeric(replicates),
                             poolSize = rep(1e+7, ncol(freq)),
                             mincov = 1,
                             MeanStart = TRUE,
                             IntGen = TRUE,
                             TA = FALSE,
                             order = 0,
                             correct = FALSE,
                             RetVal = 0)
  
  # Adjust the p-values to account for multiple comparisons
  p_vals <- p.adjust(p_vals, method = "fdr")
  
  # Returns the p-values, which is what we want from this function. P-values are not transformed to -log10
  return(p_vals)
}

getAbsPos <- function(data) {
  # Saccharomyces cerevisiae chromosome sizes based on latest release
  chromosomes <- c("C01" = 230218, #https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/
                   "C02" = 813184,
                   "C03" = 316620,
                   "C04" = 1531933,
                   "C05" = 576874,
                   "C06" = 270161,
                   "C07" = 1090940,
                   "C08" = 562643,
                   "C09" = 439888,
                   "C10" = 745751,
                   "C11" = 666816,
                   "C12" = 1078177,
                   "C13" = 924431,
                   "C14" = 784333,
                   "C15" = 1091291,
                   "C16" = 948066)
  
  # Calculates the absolute position of SNPs based on the chromosome sizes
  # For this to work, data has to be ordered by the chromosomes (i.e. all SNPs for C01, then all SNPs for C02, etc)
  ABS_POS <- c(
    data$pos[data$chr == "C01"],
    data$pos[data$chr == "C02"] + sum(chromosomes[1]),
    data$pos[data$chr == "C03"] + sum(chromosomes[1:2]),
    data$pos[data$chr == "C04"] + sum(chromosomes[1:3]),
    data$pos[data$chr == "C05"] + sum(chromosomes[1:4]),
    data$pos[data$chr == "C06"] + sum(chromosomes[1:5]),
    data$pos[data$chr == "C07"] + sum(chromosomes[1:6]),
    data$pos[data$chr == "C08"] + sum(chromosomes[1:7]),
    data$pos[data$chr == "C09"] + sum(chromosomes[1:8]),
    data$pos[data$chr == "C10"] + sum(chromosomes[1:9]),
    data$pos[data$chr == "C11"] + sum(chromosomes[1:10]),
    data$pos[data$chr == "C12"] + sum(chromosomes[1:11]),
    data$pos[data$chr == "C13"] + sum(chromosomes[1:12]),
    data$pos[data$chr == "C14"] + sum(chromosomes[1:13]),
    data$pos[data$chr == "C15"] + sum(chromosomes[1:14]),
    data$pos[data$chr == "C16"] + sum(chromosomes[1:15])
  )
  return(ABS_POS)
}

getCmhPlot <- function(data, replicates, title) {
  # Filters data to consider only the columns that have data for samples under the replicates of interest (control, low, or high)
  # The code works as it collapses rep + number + |, inducing grep to choose columns with repX OR repY or repZ
  data_filtered <- data[,grep(paste("rep", replicates, sep = "", collapse = "|"), colnames(data))]
  
  # Uses the previous function to get p-values
  data_filtered$p_vals <- getAdaptedCMH(data = data_filtered, replicates = replicates)
  
  # Creates a new column called abs_pos using the previous function
  data_filtered$abs_pos <- getAbsPos(data)
  
  # Adds the chr column - this is identical to the original dataset column, as there is no snps filtered out, only columns
  data_filtered$chr <- data$chr
  
  # Set the axis to be centered in the average position for the chromosome (average not equals median, but as we have a lot of data points, its close enough)
  axis_set <- data_filtered %>%
    group_by(chr) %>%
    summarize(center = mean(abs_pos))
  
  # Creates the plot - X axis is absolute position, colored by chromosome; y axis is -log10(p-values)
  plot(data_filtered$abs_pos,
       -log10(data_filtered$p_vals),
       type = "p",
       xaxt = "n",
       ylim = c(0,70),
       xlab = "Genomic position",
       ylab = "-log10(p-value)",
       main = title,
       col = factor(data_filtered$chr)
  )
  axis(1, at=axis_set$center, labels=axis_set$chr)
  abline(a = -log10(0.005), b = 0)
  
  # Save plot to print later
  p <- recordPlot()
  return(p)
}

# Import data from TXT
data <- read.table("data/burkephillips_exercise/SNPtable_raw.txt", header = TRUE)
data[,5:12] <- NULL # These columns are not important for this exercise

# Create vectors of different replicates
#control_replicates <- c(01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
control_replicates <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
moderate_replicates <- c(25, 26, 27, 28, 29, 30, 32, 33, 35, 36, 37, 38, 41, 42, 43, 44, 45, 46, 47, 48)
high_replicates <- c(49, 50, 51, 52, 53, 54, 55, 56, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72)

# Plots
getCmhPlot(data, replicates = control_replicates, title = "A. Control")
getCmhPlot(data, replicates = moderate_replicates, title = "B. Moderate Ethanol Treatment")
getCmhPlot(data, replicates = high_replicates, title = "C. High Ethanol Treatment")



# Attempt to concatenate everything in one big figure - FAIL!
png("burkephillips-cmh.png", width=(1200), height = (3*986))

par(mfrow = c(3, 1))
p1 <- getCmhPlot(data, replicates = control_replicates, title = "A. Control")
p2 <- getCmhPlot(data, replicates = moderate_replicates, title = "B. Moderate Ethanol Treatment")
p3 <- getCmhPlot(data, replicates = high_replicates, title = "C. High Ethanol Treatment")

grid.arrange(p1,p2,p3, nrow=3)

dev.off()

getCmhPlot(data, replicates = moderate_replicates, title = "B. Moderate Ethanol Treatment")














# Test grounds



# Ne_wp2 <- c()
# for (i in seq(1, ncol(cov), by=3)) {
#   unlisted_ne <- estimateNe(p0 = freq[,i],
#                             pt = freq[,(i+2)],
#                             cov0 = cov[,i],
#                             covt = cov[,(i+2)],
#                             t = 200,
#                             method=c("W.planII"),
#                             Ncensus=NA,
#                             poolSize=c(1e+7, 1e+7))
#   #print(as.integer(unname(unlisted_ne)))
#   
#   Ne_wp2 <- c(Ne_wp2, as.integer(unname(unlisted_ne)))
# }
# 
# Ne_pp2 <- c()
# for (i in seq(1, ncol(cov), by=3)) {
#   unlisted_ne <- estimateNe(p0 = freq[,i],
#                             pt = freq[,(i+2)],
#                             cov0 = cov[,i],
#                             covt = cov[,(i+2)],
#                             t = 200,
#                             method=c("P.planII"),
#                             Ncensus=NA,
#                             poolSize=c(1e+7, 1e+7))
#   #print(as.integer(unname(unlisted_ne)))
#   
#   Ne_pp2 <- c(Ne_pp2, as.integer(unname(unlisted_ne)))
# }
# 
# Ne_jrp2 <- c()
# for (i in seq(1, ncol(cov), by=3)) {
#   unlisted_ne <- estimateNe(p0 = freq[,i],
#                             pt = freq[,(i+2)],
#                             cov0 = cov[,i],
#                             covt = cov[,(i+2)],
#                             t = 200,
#                             method=c("JR.planII"),
#                             Ncensus=NA,
#                             poolSize=c(1e+7, 1e+7))
#   #print(as.integer(unname(unlisted_ne)))
#   
#   Ne_jrp2 <- c(Ne_jrp2, as.integer(unname(unlisted_ne)))
# }











#library(ggplot2)
#library(ggtext)



# manhplot <- ggplot(df, aes(x = abs_pos, y = -log10(p_vals), color = as.factor(chr))) +
#   geom_point(alpha = 0.75) +
#   scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
#   ylim(0,max(-log10(df$p_vals)) + 25) +
#   scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
#   scale_size_continuous(range = c(0.5,3)) +
#   labs(x = "Chromosome", y = "-log10(p)") +
#   theme_set(theme_minimal()) +
#   theme(legend.position = "none",
#         #panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   ggtitle(title)
# 
# return(manhplot)