library(ggplot2)
library(gridExtra)

# 1 - Manhattan Plot ===========================================================
GetManhattanPlot <-
  function(my_dataframe,
           Y,
           permutation_pvals = NULL,
           percentage_significance = NULL,
           title = "Manhattan Plot",
           x_label = TRUE,
           y_label = "-log10(p)",
           palette = "blue",
           y_limit_up = 100,
           y_limit_down = 0) {
  
  # This huge function plots a good manhattan plot considering all the settings I have on the arguments.
  # Some of them have been tweaked over the time to reach what I think is an optimal graph for my use case
  axis_set <-
    c("2L" = (max(my_dataframe[my_dataframe$CHROM == "2L", ]$ABS_POS) + 
                min(my_dataframe[my_dataframe$CHROM == "2L", ]$ABS_POS))/2,
      "2R" = (max(my_dataframe[my_dataframe$CHROM == "2R", ]$ABS_POS) + 
                min(my_dataframe[my_dataframe$CHROM == "2R", ]$ABS_POS))/2,
      "3L" = (max(my_dataframe[my_dataframe$CHROM == "3L", ]$ABS_POS) + 
                min(my_dataframe[my_dataframe$CHROM == "3L", ]$ABS_POS))/2,
      "3R" = (max(my_dataframe[my_dataframe$CHROM == "3R", ]$ABS_POS) + 
                min(my_dataframe[my_dataframe$CHROM == "3R", ]$ABS_POS))/2,
      "X"  = (max(my_dataframe[my_dataframe$CHROM == "X",  ]$ABS_POS) + 
                min(my_dataframe[my_dataframe$CHROM == "X",  ]$ABS_POS))/2)
  
  p <- ggplot(my_dataframe,
              aes(x = ABS_POS, 
                  y = Y,
                  color = as.factor(CHROM))) +
    geom_point(alpha = 0.75) +
    scale_x_continuous(labels = names(axis_set), breaks = axis_set) +
    {if (!x_label) scale_x_continuous(labels = NULL, breaks = axis_set)} +
    ylim(y_limit_down, y_limit_up) +
    {if (palette == "blue")
      scale_color_manual(values = rep(c("steelblue1", "steelblue4"), unique(length(axis_set))))} +
    {if (palette == "red")
      scale_color_manual(values = rep(c("firebrick1", "firebrick4"), unique(length(axis_set))))} +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = y_label) +
    theme_minimal() +
    theme(legend.position="none") +
    ggtitle(title)
  
  if (!is.null(permutation_pvals)) {
    qt90 <- quantile(as.numeric(permutation_pvals), probs = 0.90)
    qt95 <- quantile(as.numeric(permutation_pvals), probs = 0.95)
    qt99 <- quantile(as.numeric(permutation_pvals), probs = 0.99)
    qt999 <- quantile(as.numeric(permutation_pvals), probs = 0.999)
    
    p <- p +
      # 90% line
      geom_hline(yintercept = qt90, col = "darkgreen") +
      annotate("text", y=(qt90-2), x=max(my_dataframe$ABS_POS), label = "10%", col = "darkgreen") +
      # 95% line
      geom_hline(yintercept = qt95, col = "darkorange") +
      annotate("text", y=(qt95-2), x=max(my_dataframe$ABS_POS), label = "5%", col = "darkorange") +
      # 99% line
      geom_hline(yintercept = qt99, col = "burlywood4") +
      annotate("text", y=(qt99-2), x=max(my_dataframe$ABS_POS), label = "1%", col = "burlywood4") +
      # 99.9% line
      geom_hline(yintercept = qt999, col = "darkorchid") +
      annotate("text", y=(qt999-2), x=max(my_dataframe$ABS_POS), label = "0.1%", col = "darkorchid")
  }
  
  if (percentage_significance) {
    threshold <- 100 # or -log10(1e-100)
    p <- p +
      geom_hline(yintercept = threshold, col = "red")
      #annotate("text", y=(threshold-5), x=max(my_dataframe$ABS_POS), label = paste(percentage_significance*100, "% top SNPs"), col = "red")
  }
  
  return(p)
}

# 2 - Coverage =================================================================
PlotCoverage <- function(snp_table,
                         treatment,
                         gen,
                         mode = "normal") {
  
  # Mode can be "normal", "scaled", or "standardized"
  
  # Scaling function
  scaleCol <- function(val) {
    scaledCol <- val / mean(val)
    return(scaledCol)
  }
  
  # Standardize function
  standardize = function(x){
    z <- (x - mean(x)) / sd(x)
    return( z)
  }
  
  cov <- snp_table[grep(paste("N_", treatment, "_rep.._gen", gen, sep=""), colnames(snp_table))]
  
  if (mode == "scaled") {
    cov2 <- as.data.frame(apply(cov, 2, scaleCol))
  }
  
  if (mode == "standardized") {
    cov3 <- as.data.frame(apply(cov, 2, standardize))
  }
  
  plots <- list()
  for(i in 1:ncol(cov)){
    plots[[i]] <- manhplot(snp_table, Y = cov[,i], title = paste(names(cov[i]), "_", mode, sep = ""))
  }
  do.call(grid.arrange, c(plots, ncol=3, nrow=2))
}

# 3 - Allelic Trajectory =======================================================
GetAllelicTrajectoryPlot <- function(snp_table,
                                     snp,
                                     treatments = c("OB", "OBO", "nB", "nBO"),
                                     cex = 2.5) {
  
  # Define samples and filter freq dataset to contain only the interesting SNP
  snp_chrom <- gsub(":", "", unlist(strsplit(snp, ":"))[1])
  snp_pos <- as.numeric(unlist(strsplit(snp, ":"))[2])
  
  #samples <- unique(gsub("alt_", "", gsub("N_", "", colnames(data[grep("alt_|N_", colnames(data),)]))))
  freq <- snp_table[,grep("^alt_", colnames(snp_table))]/snp_table[,grep("^N_", colnames(snp_table))]
  freq$POS <- snp_table$POS
  freq$CHROM <- snp_table$CHROM
  freq_filtered <- freq[snp_table$POS == snp_pos & snp_table$CHROM == snp_chrom, ]

  
  # Transpose the dataframe to make manipulation easier
  freq_filtered_t <- as.data.frame(t(freq_filtered))
  rownames(freq_filtered_t) <- gsub("alt_", "", rownames(freq_filtered_t))
  
  # Set the axis based on what treatments I am interested in
  if (any(c("CBO", "CB", "nBO", "nB") %in% treatments)) {axis_set = c("01" = 1, "56" = 20)}
  if (any(c("EBO", "EB", "OBO", "OB") %in% treatments)) {axis_set = c("01" = 1, "20" = 20)}
  if (any(c("CBO", "CB", "nBO", "nB") %in% treatments) & any(c("EBO", "EB", "OBO", "OB") %in% treatments)) {axis_set = c("01" = 1, "20/56" = 20)}
  
  # Set the title
  title = paste("SNP ", freq_filtered_t["CHROM", ], ":", format(as.numeric(freq_filtered_t["POS", ]), big.mark = ","), sep = "")
  
  colors <- c("EB" = "blue",
              "OB" = "blue",
              "EBO" = "darkblue",
              "OBO" = "darkblue",
              "CB" = "red",
              "nB" = "red",
              "CBO" = "darkred",
              "nBO" = "darkred")
  
  strokes <- c("EB" = 2,
               "OB" = 2,
               "EBO" = 1,
               "OBO" = 1,
               "CB" = 2,
               "nB" = 2,
               "CBO" = 1,
               "nBO" = 1)
  
  # Initialize a clean plot with defined lims
  plot(1,
       type = "n",
       xlab = "", 
       ylab = "",
       xlim = c(1, 20),  
       ylim = c(0, 1),
       xaxt = "n",
       main = title,
       cex.axis = cex
  ) 
  
  # Grabs the frequency for all desired treatments
  for (treatment in treatments) {
    plotting_df <- cbind(
      freq_filtered_t[grep(paste(treatment, "_rep0._gen01", sep=""), rownames(freq_filtered_t)), 1],
      freq_filtered_t[grep(paste(treatment, "_rep.._gen20|", treatment, "_rep.._gen56", sep=""), rownames(freq_filtered_t)), 1]
    )
    
      # Gets color and stroke for treatments
    color = colors[treatment]
    stroke = strokes[treatment]
    
    # Adds a line for each treatment
    for (line in 1:nrow(plotting_df)) {
      lines(c(1,20), plotting_df[line, ], type = "b", col=color, lty=stroke)
    }
    axis(1, at=axis_set, labels=names(axis_set), cex.axis = cex)
  }
  
  # Gets the legend. Commented for now as I don't require it
  # legend(x = "center",
  #        legend = treatments,
  #        lty = strokes[treatments],
  #        col = colors[treatments],
  #        cex = cex/2)
  
  plot <- recordPlot()
  
  return(plot)
}
