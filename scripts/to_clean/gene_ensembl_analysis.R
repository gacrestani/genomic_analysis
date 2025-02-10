# Install biomaRt if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library(biomaRt)

# Get dataset
ensembl = useEnsembl(biomart="ensembl")
listDatasets(ensembl)

# Connect to Ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl")

# Define coordinates
chromosome <- "2L"
start <- 2901919
end <- 2901919

# Query for genes overlapping the coordinate
result <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
  filters = c("chromosome_name", "start", "end"),
  values = list(chromosome, start, end),
  mart = ensembl
)

print(result)
result$external_gene_name

gene_list <- c()
window <- 10000

for (i in 1:length(peaks)) {
  chromosome <- data$CHROM[data$coordinate == peaks[i]]
  pos <- data$POS[data$coordinate == peaks[i]]
  
  result <- getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
    filters = c("chromosome_name", "start", "end"),
    values = list(chromosome, pos-window, pos+window),
    mart = ensembl
  )
  
  print(result)
  gene_list <- c(gene_list, result$external_gene_name)
}

#gene_list_backup <- gene_list
gene_list <- gene_list_backup

filtered_vector <- my_vector[!grepl(pattern, my_vector)]

gene_list <- gene_list[!grepl("^lncRNA*", gene_list)]
gene_list_unique <- unique(gene_list)

cat(gene_list_unique, sep = "\n")

