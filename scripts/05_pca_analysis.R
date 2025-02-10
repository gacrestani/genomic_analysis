library(ggplot2)
library(ggpubr)
library(ggrepel)

PreparePca <- function(snp_table) {
  
  freq <- GetFreq(snp_table)
  pca_df <- t(as.matrix(freq)) # Transpose the data frame to fit prcomp standard
  pca <- prcomp(pca_df)
  
  pca_data <- data.frame(sample=gsub("alt_|N_", "", colnames(freq)),
                         X = pca$x[,1],
                         Y = pca$x[,2],
                         Z = pca$x[,3]) # Select first three Principal Components
  
  # Create grouping variable
  pca_data$population <- factor(
    gsub("([A-Z]+)_rep.._(gen..)",
         "\\1_\\2",
         pca_data$sample))
  
  pca_data$variance <- pca$sdev^2 # Create var column
  pca_data$variance_percentage <-
    round(pca_data$variance / sum(pca_data$variance)*100, 2)
  
  clustering_result <- kmeans(pca_df, centers = 3, nstart = 25)
  pca_data$cluster <- as.factor(clustering_result$cluster)
  
  return(pca_data)
}

PlotPca <- function(pca_data,
                    label = TRUE,
                    title = NULL) {
  
  # If label = TRUE, all data points will be labeled independently
  # If label = FALSE, data points will be colored by population
  
  # Pick a palette that fits your data
  rbpalette=c("red",
              "darkred",
              "magenta",
              "darkorchid4",
              "blue",
              "darkblue",
              "green",
              "darkgreen")
  
  # This one has the text label for all the samples
  pca_plot <-
    ggplot(data = pca_data,
           aes(x=X,
               y=Y,
               label=sample,
               color=population,
               shape = cluster)) +
    geom_point(size = 3) +
    scale_color_manual(values = rbpalette) +
    xlab(paste("PC1 - ", pca_data$var.per[1], "%", sep="")) +
    ylab(paste("PC2 - ", pca_data$var.per[2], "%", sep="")) +
    ggtitle(title)  +
    {if (label) geom_text_repel()} +
    {if (label) theme(legend.position="none")} +
    theme_bw()
  
  return(pca_plot)
}