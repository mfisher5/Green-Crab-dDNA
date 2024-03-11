eDNA_index <- function(x) { #where x is a dataframe with taxa/OTUs/etc in rows, and samples in columns
  rowMax <- function(x){apply(x, MARGIN = 1, FUN = max)}
  temp <- sweep(x, MARGIN = 2, STATS = colSums(x), FUN = "/") #create proportion for each taxon/OTU within each sample
  sweep(temp, MARGIN = 1, STATS = rowMax(temp), FUN = "/")
}
