filter.BEDPE <- function(bed, filter, redundant = F){
  require(GenomicRanges)
  colnames(filter) <- c('chr', 'start', 'end')
  filter <- makeGRangesFromDataFrame(filter)
  
  bed1 <- bed[, 1:3]; colnames(bed1) <- c('chr', 'start', 'end')
  bed1 <- makeGRangesFromDataFrame(bed1)
  ov1 <- findOverlaps(bed1, filter)
  
  if(redundant){
    finders <- queryHits(ov1)
  }else{
    bed2 <- bed[, 4:6]; colnames(bed2) <- c('chr', 'start', 'end')
    bed2 <- makeGRangesFromDataFrame(bed1)
    ov2 <- findOverlaps(bed2, filter)
    
    finders <- unique(c(queryHits(ov1),queryHits(ov2)))
  }
  
  keepers <- bed[-finders, ]
  return(keepers)
}
