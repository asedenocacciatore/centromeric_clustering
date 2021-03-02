#' Exclude some regions from a Hi-C object (by removing that info)
#'
#' @param exp GENOVA (<v1.0) data structure for the experiment of interest
#' @param chromsToUse a vector containing the chromosomes that should be considered
#' @param excludeBed regions to be excluded in BED format
#'
#' @return dataframe similar to ICE without information for the excluded bins
#' 
filter.ICE <- function(exp, chromsToUse, excludeBed = NULL){
  
  # Subset ICE data to only consider interactions between not-excluded regions
  finders <- exp$ABS$V4[exp$ABS$V1 %in% chromsToUse]
  
  if(!is.null(excludeBed)){
    exclude <- c()
    for(row in 1:nrow(excludeBed)){
      # This can probably done more elegantly, but apply is not working as expected
      exclude <- c(exclude, exp$ABS$V4[exp$ABS$V1 == excludeBed[row, 1] & exp$ABS$V2 >= excludeBed[row, 2] & exp$ABS$V3 <= excludeBed[row, 3]])
    }
    finders <- finders[!finders %in% exclude] 
  }
  
  keepers <- exp$ICE$V1 %in% finders & exp$ICE$V2 %in% finders
  dat <- exp$ICE[keepers, ]
  return(dat)
}
