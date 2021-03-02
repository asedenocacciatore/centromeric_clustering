MeasureDistanceCentromere <- function(gr, centromeres) {
  # For a GRanges object and centromere information, calculate distance to centromere
  # Obtained a scaled distance in relation to chromosome arm length
  
  # Calculate distance to centromere
  gr$distance_to_centromere <- mcols(distanceToNearest(gr, centromeres))$distance
  gr$proportion_to_centromere <- as.numeric(gr$distance_to_centromere)
  
  # For each chromosome
  #   p-arm: normalize the distance by the first base - start centromere distance
  #   q-arm: normalize by the distance by the end centromere - last base
  for (chr in seqlevels(gr)) {
    centromeres.chr <- centromeres[seqnames(centromeres) == chr]
    
    left_arm <- start(centromeres.chr)
    right_arm <- max(end(gr[seqnames(gr) == chr, ]) - end(centromeres.chr))
    
    idx <- seqnames(gr) == chr
    gr.tmp <- gr[idx]
    
    norm_dis <- ifelse(start(gr.tmp) < start(centromeres.chr),
                       gr.tmp$distance_to_centromere / left_arm,
                       gr.tmp$distance_to_centromere / right_arm)

    mcols(gr[idx])[, 'proportion_to_centromere'] <- norm_dis
  }
  
  # Return gr with new info
  return(gr)
}
