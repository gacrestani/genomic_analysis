source("functions.R")

df <- importData()

snpeff_df <- fread(file = "/nfs3/IB/Burke_Lab/Crestani/nextflow/results/annotated_snps.txt", header = TRUE)
snpeff_df <- addAbsPos(snpeff_df)
snpeff_df <- semi_join(snpeff_df, df, by = "ABS_POS")

cmh.pvals <- fread("/nfs3/IB/Burke_Lab/Crestani/flylong_project/genomic_analyses/results/cmh_pvals.csv")
perm.pvals <- fread("/nfs3/IB/Burke_Lab/Crestani/flylong_project/genomic_analyses/results/perm_pvals.csv")

df$o_pvals <- cmh.pvals$o
df$ANN <- snpeff_df$ANN

qt90 <- quantile(as.numeric(perm.pvals$o), probs = 0.90)
qt95 <- quantile(as.numeric(perm.pvals$o), probs = 0.95)
qt99 <- quantile(as.numeric(perm.pvals$o), probs = 0.99)
qt999 <- quantile(as.numeric(perm.pvals$o), probs = 0.999)

###

cutoff <- qt95

df_filtered <- df[df$o_pvals > cutoff,]

vector <- c()

for (i in 1:nrow(df_filtered)) {
  #print(i)
  a <- strsplit(as.character(df_filtered$ANN), ",")[[i]]
  #print(a[1])
  vector <- append(vector, a[1])
  
}

ann_df <- as.data.frame(do.call(rbind, strsplit(as.character(vector), "\\|")))

cat(unique(ann_df$V5[ann_df$V8 == "protein_coding"]), sep = "\n")
cat(unique(ann_df$V5), sep = "\n")


