## Interchromosomal contacts enrichment
chromosomeMatrix.plot <- function(mat, chromsToOmit = NULL){
  # mat = result from chromosomeMatrix
  if(!is.null(chromsToOmit)){
    keep <- !colnames(mat$normMat) %in% chromsToOmit
    keep <- outer(keep, keep, '*')
    mat$rawCounts <- mat$rawCounts*keep
    mat$normMat <- mat$normMat*keep
  }
  
  totalInter <- sum(mat$rawCounts[upper.tri(mat$rawCounts)])
  chromInter <- (colSums(mat$rawCounts) - diag(mat$rawCounts))/totalInter
  mat$expected <- outer(chromInter, chromInter, '*')*totalInter
  
  oe <- mat$rawCounts/mat$expected
  rownames(oe) <- rownames(mat$normMat); colnames(oe) <- colnames(mat$normMat)
  diag(oe) <- NA
  long.oe <- melt(oe)
  ggplot(long.oe, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_x_discrete(labels = function(x) format(sub('chr', '', x))) +
    scale_y_discrete(limits = rev(levels(long.oe$Var1)), labels = function(x) format(sub('chr', '', x))) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1, space = "Lab", name="O/E") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_fixed() +
    labs(x="Human chromosomes", y="Human chromosomes")
}

# Hi-C trans maps of translocated chromosomes
select.trans <- function(exp, chrom1, chrom2, na.rm=F){
  #id1 <- exp$ABS[exp$ICE$V1, 1]
  #id2 <- exp$ABS[exp$ICE$V2, 1]
  #sel <- id1 %in% c(chrom1, chrom2) & id2 %in% c(chrom1, chrom2) & id1 != id2
  
  idx1 <- exp$ABS[exp$ABS$V1 == chrom1, 4]
  idx2 <- exp$ABS[exp$ABS$V1 == chrom2, 4]
  
  x <- rep(idx1, length(idx2))
  y <- rep(idx2, each=length(idx1))
  mat <- exp$ICE[list(x, y)]
  if(na.rm){
    mat <- mat[!is.na(mat$V3), ]
  }

  return(mat)
}

hic_transplot <- function(exp, chrom1, chrom2){
  # Order chroms to compare
  chroms <- unique(exp$ABS$V1)
  if(which(chroms == chrom1) > which(chroms == chrom2)){
    tmp <- chrom1
    chrom1 <- chrom2
    chrom2 <- tmp
    rm(tmp)
  }
  
  # Select the contacts between the two selected chroms
  mat <- select.trans(exp, chrom1, chrom2)
  
  z <- c(quantile(na.exclude(mat$V3), 0.00),
         quantile(na.exclude(mat$V3), 0.99))
  th <- quantile(na.omit(mat$V3), 0.99)
  mat[mat$V3 > th, 3] <- th
  
  # Plot
  ggplot(mat, aes(x=exp$ABS[V1, 2], y=exp$ABS[V2, 2], fill=V3)) +
    geom_raster() +
    scale_fill_gradientn(colours = c('white', 'red', '#990000'), limits = z, name='Contacts', na.value = 'white') +
    scale_x_discrete(labels=function(x) format(x/1e6),
                     limits=c(0, max(exp$ABS[mat$V1, 2]))) +
    scale_y_discrete(labels=function(x) format(x/1e6),
                     limits=c(0, max(exp$ABS[mat$V2, 2]))) +
    labs(x=gsub('chr', 'Chr ', chrom1), y=gsub('chr', 'Chr ', chrom2)) +
    coord_fixed() +
    theme_bw()
}
