# Sliding window
getmeans <- function(x,mafs,pos){
  half.window.size <<-50000   ### you can play with this number, in this case the window size is 100kb
  window <- (pos > (x - half.window.size)) & (pos < (x + half.window.size))
  N <<- sum(window)
  Y <<- as.matrix(mafs[window],ncol=1)
  out <- mean(Y)
}

# So, imagine you're working with chrX for the B1 population (for the sample code below, I'm assuming the subsetted SNP table has the name "chrX" but obviously this can be adjusted). You'd want to define a vector of positions for the windows, then apply the function.
min=min(data$pos)
max=max(data$pos)
testpositions=seq(min+100000,max-100000,50000)

# this vector of positions specifies the step size (in this case, 100kb to impose non-overlap, but use a smaller number e.g. 2000 if you want overlapping windows)

sw_B_Fst <- sapply(testpositions,function(x) getmeans(x,df$F_st,df$pos))
sw_B_Het <- sapply(testpositions,function(x) getmeans(x,df$H_t,df$pos))

# this assumes that the B1het is a column in the chrX data frame with the B1 heterozygosity values


#Let me know if you have questions or if something doesn't look right!  Also, when you are plotting, I recommend dividing positions by 10^6 - this converts base pairs to megabases (MB).  This will look nicer when you plot.  So for example:

KB <- testpositions/10^3

plot(KB,sw_B_Fst)
plot(KB,sw_B_Het)


## Another version
# Sliding window settings ####
getMeans <- function(x,mafs,pos){
  half.window.size <<- 50000 # you can play with this number, in this case the window size is 50kb
  window <- (pos > (x - half.window.size)) & (pos < (x + half.window.size))
  N <<- sum(window)
  Y <<- as.matrix(mafs[window],ncol=1)
  out <- mean(Y)
}

min=min(data$POS)
max=max(data$POS)
testpositions=seq(min+100000,max-100000,50000)
KB <- testpositions/10^3

p_vals <- getAdaptedCMH(data = data, freq_df = freq_df, treatment = "CBO", log = TRUE)
p_vals_smooth <- sapply(testpositions, function(x) getMeans(x, p_vals, data$ABS_POS))
plot(KB, p_vals_smooth, type = "l", main = "CBO")