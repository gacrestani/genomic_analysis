# Data Import #####
require(tidyverse)
require(ggplot2)
require(reshape2)

# Reads SNPs data
data <- read.table("data/graves/SNP_table.txt", header = TRUE)

# Filters for only X chromossome
data <- data[data$chr == "X", ]

# Define treatments
treatments <- c("NCO", "AO", "BO", "ACO","B", "CO")

# Filters for only regular chromossomes
#data <- data[data$chr == c("X","2l", "2R", "3L", "3R", "4"), ]
#####



# Frequencies and Fst #####
getFrequencies <- function(replicate, treatment, data = data) {
  # Iterates over treatments and replicates to calculate p, q and heterozygosity
  iteration_colname <- paste("rep", replicate, "_", treatment, "_", sep = "")
  df <- data.frame(matrix(ncol=0, nrow=nrow(data)))
  df[paste(iteration_colname, "p", sep = "")] <- data[,paste(iteration_colname, "1", sep = "")] / (data[,paste(iteration_colname, "1", sep = "")] + data[,paste(iteration_colname, "0", sep = "")])
  df[paste(iteration_colname, "q", sep = "")] <- 1 - df[paste(iteration_colname, "p", sep = "")]
  df[paste(iteration_colname, "2pq", sep = "")] <- 2 * df[paste(iteration_colname, "p", sep = "")] * df[paste(iteration_colname, "q", sep = "")]
  return(df)
}

getFst <- function(replicate, treatment, data) {
  treatment_colname <- paste()
  iteration_colname <- paste("rep", replicate, "_", treatment, "_", sep = "")
  
  df <- data.frame(matrix(ncol=0, nrow=nrow(data)))
  
  df[paste(treatment, "_p_t", sep = "")] <- 0
  for (i in 1:replicate) {
    df[paste(treatment, "_p_t", sep = "")] <- df[paste(treatment, "_p_t", sep = "")] + data[paste(iteration_colname, "p", sep = "")]
  }
  df[paste(treatment, "_p_t", sep = "")] <- df[paste(treatment, "_p_t", sep = "")]/replicate
  
  df[paste(treatment, "_q_t", sep = "")] <- 1 - df[paste(treatment, "_p_t", sep = "")]
  df[paste(treatment, "_2pq_t", sep = "")] <- 2 * df[paste(treatment, "_p_t", sep = "")] * df[paste(treatment, "_q_t", sep = "")]
  
  df[paste(treatment, "_Hs", sep = "")] <- 0
  for (i in 1:replicate) {
    df[paste(treatment, "_Hs", sep = "")] <- df[paste(treatment, "_Hs", sep = "")] + data[paste(iteration_colname, "2pq", sep = "")]
  }
  df[paste(treatment, "_Hs", sep = "")] <- df[paste(treatment, "_p_t", sep = "")]/replicate
  
  df[paste(treatment, "Fst", sep = "")] <-  df[paste(treatment, "_2pq_t", sep = "")] - df[paste(treatment, "_Hs", sep = "")]/df[paste(treatment, "_2pq_t", sep = "")]
  
  return(df)
}

# Create frequencies df and Fst df
freq_df <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
for (i in 1:5) {
  freq_df <- cbind(freq_df, getFrequencies(i, treatment = treatments, data))
}

Fst_df <- getFst(replicate = 5, treatment = treatments, data = freq_df) # This has some -inf

# Plot raw heterozigosity
par(mfrow=c(5,1))


raw_hetero <- freq_df[,grep("_2pq", colnames(freq_df))]
raw_hetero$pos <- data$pos

het_melt <- melt(raw_hetero, id=c("pos"))

ggplot() +
  geom_line(data = het_melt %>% filter(str_detect(variable, "NCO")), aes(x=pos, y=value, group=factor(variable)))
#####



# PCA for p #####
pca_df <- freq_df[,grep("_p", colnames(freq_df))]
pca_df <- cbind(data[1:2], pca_df)

pca <- prcomp(t(pca_df[,3:ncol(pca_df)]), scale = TRUE)

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 2)

pca.data <- data.frame(Sample=rownames(pca$x), X = pca$x[,1], Y = pca$x[,2])

pca.data$pop <- substring(pca.data$Sample, 6, 7)

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, color = as.factor(pca.data$pop))) +
  geom_point() + 
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))
#####



#### Sliding window
getmeans <- function(x,mafs,pos){
  half.window.size <<- 50000   ### you can play with this number, in this case the window size is 100kb
  window <- (pos > (x - half.window.size)) & (pos < (x + half.window.size))
  N <<- sum(window)
  Y <<- as.matrix(mafs[window],ncol=1)
  out <- mean(Y)
}

# So, imagine you're working with chrX for the B1 population (for the sample code below, I'm assuming the subsetted SNP table has the name "chrX" but obviously this can be adjusted).  You'd want to define a vector of positions for the windows, then apply the function.
min=min(data$pos)
max=max(data$pos)
testpositions=seq(min+100000,max-100000,100000)

# this vector of positions specifies the step size (in this case, 100kb to impose non-overlap, but use a smaller number e.g. 2000 if you want overlapping windows)
df$F_st[is.nan(df$F_st)] <- 0

sw_B_Fst <- sapply(testpositions,function(x) getmeans(x,df$F_st,df$pos))
sw_B_Het <- sapply(testpositions,function(x) getmeans(x,df$H_t,df$pos))

# this assumes that the B1het is a column in the chrX data frame with the B1 heterozygosity values


#Let me know if you have questions or if something doesn't look right!  Also, when you are plotting, I recommend dividing positions by 10^6 - this converts base pairs to megabases (MB).  This will look nicer when you plot.  So for example:

KB <- testpositions/10^3

# PLOTS
#plot(KB,sw_B_Fst, type = "l")
#plot(KB,sw_B_Het, type="l")
####


### Perform Fisher's exact and store the p-values

# Comparisons
## ACO1 vs. CO1
## ACO1 vs. B1
## CO1 vs. B1

getFisherP <- function(samples, main) {
  filtered_snps <- data[samples]
  a<-as.matrix(filtered_snps[, c(1:2)])
  b<-as.matrix(filtered_snps[, c(3:4)])
  
  p_vals <- c()
  for (i in 1:nrow(data)) {
    result <- fisher.test(rbind(a[i,], b[i,]))
    p_vals <- c(p_vals, result$p.value)
  }
  
  return(p_vals)
}

acoco <- getFisherP(c("rep1_ACO_1", "rep1_ACO_0", "rep1_CO_1", "rep1_CO_0"))
acob <- getFisherP(c("rep1_ACO_1", "rep1_ACO_0", "rep1_B_1", "rep1_B_0"))
cob <- getFisherP(c("rep1_CO_1", "rep1_CO_0", "rep1_B_1", "rep1_B_0"))

acoco <- -log10(acoco)
acob <- -log10(acob)
cob <- -log10(cob)

sw_B_pval1 <- sapply(testpositions,function(x) getmeans(x,acoco,data$pos))
sw_B_pval2 <- sapply(testpositions,function(x) getmeans(x,acob,data$pos))
sw_B_pval3 <- sapply(testpositions,function(x) getmeans(x,cob,data$pos))

par(mfrow=c(3,1))
plot(KB,sw_B_pval1, type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "ACO1 vs CO1")
plot(KB,sw_B_pval2, type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "ACO1 vs B1")
plot(KB,sw_B_pval3, type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "CO1 vs B1")

# CMH
getCmhP <- function(samples, main) {
  filtered_snps <- data[samples]
  a<-as.matrix(filtered_snps[, c(1:2)])
  b<-as.matrix(filtered_snps[, c(3:4)])
  
  p_vals <- c()
  for (i in 1:nrow(data)) {
    result <- mantelhaen.test(a[i,], b[i,])
    p_vals <- c(p_vals, result$p.value)
  }
  
  return(p_vals)
}


#renv::install("ThomasTaus/poolSeq")
library(poolSeq)

## LOG10 BEFORE
# ACO vs CO
acoA <- t(data[ ,c("rep1_ACO_1", "rep2_ACO_1", "rep3_ACO_1", "rep4_ACO_1")])
acoa <- t(data[ ,c("rep1_ACO_0", "rep2_ACO_0", "rep3_ACO_0", "rep4_ACO_0")])
coA <- t(data[ ,c("rep1_CO_1", "rep2_CO_1", "rep3_CO_1", "rep4_CO_1")])
coa <- t(data[ ,c("rep1_CO_0", "rep2_CO_0", "rep3_CO_0", "rep4_CO_0")])
ACOvsCO <- cmh.test(acoA, acoa, coA, coa, log = TRUE)

# ACO vs B
acoA <- t(data[ ,c("rep1_ACO_1", "rep2_ACO_1", "rep3_ACO_1", "rep4_ACO_1")])
acoa <- t(data[ ,c("rep1_ACO_0", "rep2_ACO_0", "rep3_ACO_0", "rep4_ACO_0")])
bA <- t(data[ ,c("rep1_B_1", "rep2_B_1", "rep3_B_1", "rep4_B_1")])
ba <- t(data[ ,c("rep1_B_0", "rep2_B_0", "rep3_B_0", "rep4_B_0")])
ACOvsB <- cmh.test(acoA, acoa, bA, ba, log = TRUE)

# B vs CO
bA <- t(data[ ,c("rep1_B_1", "rep2_B_1", "rep3_B_1", "rep4_B_1")])
ba <- t(data[ ,c("rep1_B_0", "rep2_B_0", "rep3_B_0", "rep4_B_0")])
coA <- t(data[ ,c("rep1_CO_1", "rep2_CO_1", "rep3_CO_1", "rep4_CO_1")])
coa <- t(data[ ,c("rep1_CO_0", "rep2_CO_0", "rep3_CO_0", "rep4_CO_0")])
BvsCO <- cmh.test(coA, coa, bA, ba, log = TRUE)

plot(data$pos,ACOvsCO, type = "l", ylab = "-log10 pval", xlab = "")
abline(h=54.54)

# Plot
par(mfrow=c(3,1))
plot(data$pos,ACOvsCO, type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "ACO vs CO")
plot(data$pos,ACOvsB, type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "ACO vs B")
plot(data$pos,BvsCO, type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "CO vs B")



## LOG 10 AFTER
# ACO vs CO
acoA <- t(data[ ,c("rep1_ACO_1", "rep2_ACO_1", "rep3_ACO_1", "rep4_ACO_1")])
acoa <- t(data[ ,c("rep1_ACO_0", "rep2_ACO_0", "rep3_ACO_0", "rep4_ACO_0")])
coA <- t(data[ ,c("rep1_CO_1", "rep2_CO_1", "rep3_CO_1", "rep4_CO_1")])
coa <- t(data[ ,c("rep1_CO_0", "rep2_CO_0", "rep3_CO_0", "rep4_CO_0")])
ACOvsCO <- cmh.test(acoA, acoa, coA, coa)

# ACO vs B
acoA <- t(data[ ,c("rep1_ACO_1", "rep2_ACO_1", "rep3_ACO_1", "rep4_ACO_1")])
acoa <- t(data[ ,c("rep1_ACO_0", "rep2_ACO_0", "rep3_ACO_0", "rep4_ACO_0")])
bA <- t(data[ ,c("rep1_B_1", "rep2_B_1", "rep3_B_1", "rep4_B_1")])
ba <- t(data[ ,c("rep1_B_0", "rep2_B_0", "rep3_B_0", "rep4_B_0")])
ACOvsB <- cmh.test(acoA, acoa, bA, ba)

# B vs CO
bA <- t(data[ ,c("rep1_B_1", "rep2_B_1", "rep3_B_1", "rep4_B_1")])
ba <- t(data[ ,c("rep1_B_0", "rep2_B_0", "rep3_B_0", "rep4_B_0")])
coA <- t(data[ ,c("rep1_CO_1", "rep2_CO_1", "rep3_CO_1", "rep4_CO_1")])
coa <- t(data[ ,c("rep1_CO_0", "rep2_CO_0", "rep3_CO_0", "rep4_CO_0")])
BvsCO <- cmh.test(coA, coa, bA, ba)

# Plot
par(mfrow=c(3,1))
plot(data$pos,-log10(ACOvsCO), type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "ACO vs CO")
plot(data$pos,-log10(ACOvsB), type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "ACO vs B")
plot(data$pos,-log10(BvsCO), type = "l", ylab = "-log10 pval", xlab = "")
legend("topright", legend = "CO vs B")


# Sliding window
getmeans <- function(x,mafs,pos){
  half.window.size <<- 5000   ### you can play with this number, in this case the window size is 100kb
  window <- (pos > (x - half.window.size)) & (pos < (x + half.window.size))
  N <<- sum(window)
  Y <<- as.matrix(mafs[window],ncol=1)
  out <- mean(Y)
}

# So, imagine you're working with chrX for the B1 population (for the sample code below, I'm assuming the subsetted SNP table has the name "chrX" but obviously this can be adjusted).  You'd want to define a vector of positions for the windows, then apply the function.
min=min(data$pos)
max=max(data$pos)
testpositions=seq(min+10000,max-10000,2000)
KB <- testpositions/10^3

sw_ACOvsCO <- sapply(testpositions,function(x) getmeans(x,ACOvsCO,data$pos))
sw_ACOvsB <- sapply(testpositions,function(x) getmeans(x,ACOvsB,data$pos))
sw_BvsCO <- sapply(testpositions,function(x) getmeans(x,BvsCO,data$pos))

# Plot
par(mfrow=c(3,1))
plot(KB,sw_ACOvsCO, type = "l", ylab = "-log10 pval", xlab = "", col="red")
legend("topright", legend = "ACO vs CO")
plot(KB,sw_ACOvsB, type = "l", ylab = "-log10 pval", xlab = "", col="red")
legend("topright", legend = "ACO vs B")
plot(KB,sw_BvsCO, type = "l", ylab = "-log10 pval", xlab = "", col="red")
legend("topright", legend = "CO vs B")

# Plot
par(mfrow=c(3,1))
plot(KB,-log10(sw_ACOvsCO), type = "l", ylab = "-log10 pval", xlab = "", col="red")
legend("topright", legend = "ACO vs CO")
plot(KB,-log10(sw_ACOvsB), type = "l", ylab = "-log10 pval", xlab = "", col="red")
legend("topright", legend = "ACO vs B")
plot(KB,-log10(sw_BvsCO), type = "l", ylab = "-log10 pval", xlab = "", col="red")
legend("topright", legend = "CO vs B")


# Add chromossomes
# Sizes
chrX <- 19569407
chr2L <- 22940760
chr2R <- 25229534
chr3L <- 20692330
chr3R <- 29815060
chr4 <- 1113317

# Release FB2023_01 vs FB2013_03
chrX > max(data$pos[data$chr == "X"]) # False, 19,569,407 vs 22,421,823
chr2L > max(data$pos[data$chr == "2L"]) # False, 22,940,760 vs 22,982,697
chr2R > max(data$pos[data$chr == "2R"]) # True, 25,229,534 vs 21,145,059 
chr3L > max(data$pos[data$chr == "3L"]) # False, 20,692,330 vs 24,535,547
chr3R > max(data$pos[data$chr == "3R"]) # True, 29,815,060 vs 27,893,095
chr4 > max(data$pos[data$chr == "4"]) # False, 1,113,317 vs 1,275,692

dataX <- data$pos[data$chr == "X", ]
data2L <- data[data$chr == "2L", ]
data2R <- data[data$chr == "2R", ]
data3L <- data[data$chr == "3L", ]
data3R <- data[data$chr == "3R", ]
data4 <- data[data$chr == "4", ]

# Randomize replicates - ACO vs CO

iterations <- 1000

for (i in 1:iterations) {
  choices <- c("rep1_ACO", "rep2_ACO", "rep3_ACO", "rep4_ACO", "rep5_ACO", "rep1_CO", "rep2_CO", "rep3_CO", "rep4_CO", "rep5_CO")
  choices_rand <- sample(choices, length(choices))
  
  acoA <- t(data[ ,c(paste(choices_rand[1], "_1", sep = ""),
                     paste(choices_rand[2], "_1", sep = ""),
                     paste(choices_rand[3], "_1", sep = ""),
                     paste(choices_rand[4], "_1", sep = ""),
                     paste(choices_rand[5], "_1", sep = ""))])
  
  acoa <- t(data[ ,c(paste(choices_rand[1], "_0", sep = ""),
                     paste(choices_rand[2], "_0", sep = ""),
                     paste(choices_rand[3], "_0", sep = ""),
                     paste(choices_rand[4], "_0", sep = ""),
                     paste(choices_rand[5], "_0", sep = ""))])
  
  coA <- t(data[ ,c(paste(choices_rand[6], "_1", sep = ""),
                    paste(choices_rand[7], "_1", sep = ""),
                    paste(choices_rand[8], "_1", sep = ""),
                    paste(choices_rand[9], "_1", sep = ""),
                    paste(choices_rand[10], "_1", sep = ""))])
  
  coa <- t(data[ ,c(paste(choices_rand[6], "_0", sep = ""),
                    paste(choices_rand[7], "_0", sep = ""),
                    paste(choices_rand[8], "_0", sep = ""),
                    paste(choices_rand[9], "_0", sep = ""),
                    paste(choices_rand[10], "_0", sep = ""))])
  
  ACOvsCO_rand <- cmh.test(acoA, acoa, coA, coa, log = TRUE)
  ACOvsCO_rand <- na.omit(ACOvsCO_rand)
  
  p_val_list <- append(p_val_list, max(ACOvsCO_rand))  
}

p_val_list <- list()

hist(as.numeric(p_val_list))

p_val_list_iteration1 <- p_val_list
# Should we remove that position from the dataset?

#qt(0.05, as.numeric(p_val_list))

quantile(as.numeric(p_val_list), probs = 0.90)

plot(data$pos,ACOvsCO, type = "p", ylab = "-log10 pval", xlab = "")
abline(h=54.54, col="blue")
abline(h=49.89, col="red")

data <- na.omit(data)

data$pval_aco_co <- ACOvsCO
data[max(data$pval_aco_co), "pos"]
