#' eDNA index
#'
#' `eDNAINDEX` calculates the eDNA index.
#'
#' This function was written by Ryan Kelly's lab at UW 
#' 
#' @param x  a dataframe with taxa/OTUs/etc in rows, and samples in columns


eDNAINDEX <- function(x) {
  rowMax <- function(x){apply(x, MARGIN = 1, FUN = max)}
  temp <- sweep(x, MARGIN = 2, STATS = colSums(x), FUN = "/") #create proportion for each taxon/OTU within each sample
  sweep(temp, MARGIN = 1, STATS = rowMax(temp), FUN = "/")
}