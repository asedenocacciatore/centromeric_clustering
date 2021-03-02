#' Get proportion of trans interactions per bin
#'
#' @param exp Experiment object as created by GENOVA::construct.experiment
#' @param chromsToUse List of chromosomes that should be included
#' @param excludeBed Data.frame in BED format for regions that should be filtered
#' @param interArm If TRUE then instead of trans/cis contacts per bin, the proportion inter-arm/cis contacts will be computed
#'
#' @return Data.frame with sample, chromosome, bins and trans.ratio per bin
#' 
score.trans <- function(exp, chromsToUse = NULL, excludeBed = NULL, interArm = F){
  
  if(interArm) message("Computing proportion of interarm contacts") else message("Computing proportion of trans contacts") 
  
  if(is.null(chromsToUse)){
    chromsToUse <- exp$CHRS
    dat <- exp$ICE
  }else{
    dat <- filter.ICE(exp, chromsToUse, excludeBed) # consider only interactions between chromsToUse
  }
  
  trans.percentages <- data.frame()
  
  #select the chromosome names in a specific order
  exp$ABS[,1] <- as.character(exp$ABS[,1])
  
  for(chr in chromsToUse){
    message(chr)
    #select the chromosome ids from the HiC indexes
    chr.id <- exp$ABS[exp$ABS$V1 == chr, 4]
    
    # get all cis interactions of that chromosome (as vectors is faster)
    all <- dat$V1 %in% chr.id | dat$V2 %in% chr.id & dat$V1 != dat$V2 # remove diagonal to avoid counting it twice
    pos1 <- dat$V1[all]; pos2 <- dat$V2[all]; val <- dat$V3[all]
    cis <- pos1 %in% chr.id & pos2 %in% chr.id
    cis.sum <- tapply(c(val[cis], val[cis]), c(pos1[cis], pos2[cis]), sum, na.rm = T)
    
    diag <- dat$V1 %in% chr.id & dat$V1 == dat$V2
    diagonal <- dat$V3[diag]
    names(diagonal) <- dat$V1[diag]
    cis.sum <- tapply(c(cis.sum, diagonal), c(names(cis.sum), names(diagonal)), sum) # add diagonal
    
    indexes <- chr.id[!chr.id %in% names(cis.sum) & chr.id %in% unique(c(dat$V1, dat$V2))]
    cis.sum <- c(cis.sum, setNames(rep(0, length(indexes)), as.list(indexes))) # add indexes for which there are no cis contacts
    cis.sum <- cis.sum[as.character(chr.id)]
    
    if(interArm){
      # find centromeres on the cis contacts
      tmp <- cis.sum %>% is.na() %>% as.vector() %>% rle()
      Ci <- match(max(tmp$lengths[tmp$values == T]), tmp$lengths)
      start <- max(cumsum(tmp$lengths)[Ci-1], 0) + 1 # for cases where Ci == 0, it starts at the 1st bin
      end <- cumsum(tmp$lengths)[Ci]
      
      # get the inter-arm interactions
      p <- if(start > 1) chr.id[1:start] else NA
      q <- chr.id[end:length(chr.id)]
        # skip when there is no contact data for one of the arms
      if(!any(c(dat$V1, dat$V2) %in% p) | !any(q %in% c(dat$V1, dat$V2))) next
        
      inter <- (pos1 %in% p & pos2 %in% q) | (pos1 %in% q & pos2 %in% p)
      inter.sum <- tapply(c(val[inter], val[inter]), c(pos1[inter], pos2[inter]), sum, na.rm = T)
      inter.sum <- inter.sum[as.character(chr.id)]
      
      trans.percentages <- rbind(trans.percentages, data.frame(sample = exp$NAME, chrom = chr, 
                                                                     bin = 1:length(cis.sum)*exp$RES/2, 
                                                                     cis = cis.sum,
                                                                     trans= inter.sum,
                                                                     trans.ratio = inter.sum/(cis.sum)))
    }else{
      # get trans interactions
      trans.sum <- tapply(c(val[!cis], val[!cis]), c(pos1[!cis], pos2[!cis]), sum, na.rm = T)
      trans.sum <- trans.sum[as.character(chr.id)]
      
      trans.percentages <- rbind(trans.percentages, data.frame(sample = exp$NAME, chrom = chr, 
                                                               bin = 1:length(cis.sum)*exp$RES/2, 
                                                               cis = cis.sum,
                                                               trans = trans.sum,
                                                               trans.ratio = trans.sum/(cis.sum+trans.sum)))
    }
  }
  return(trans.percentages)
}

compare.transScore <- function(expList, chromsToUse, excludeBed = NULL, control='WT', interArm = F){
  require(tidyr)
  require(dplyr)
  
  trans <- data.frame()
  for(exp in expList){
    message(exp$NAME)
    exp_score <- score.trans(exp, chromsToUse, excludeBed, interArm)
    trans <- rbind(trans, exp_score)
  }
  trans$sample <- as.factor(trans$sample)
  trans_wide <- spread(trans %>% select(-cis, -trans), sample, trans.ratio)
  
  samples <- unique(trans$sample)
  samples <- match(samples, colnames(trans_wide))
  control <- match(control, colnames(trans_wide))
  trans_wide[, samples] <- log2(trans_wide[, samples]/trans_wide[, control])
  return(trans_wide)
}

plot.transScore <- function(transScore, control = 'WT', plot.type = 'pileup', combine.arms = T, interArm = F){
  require(tidyr)
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  
  if(interArm){ mode <- 'interarm' }else{ mode <- 'trans'}
  
  chroms <- unique(transScore$chrom)
  for (chr in chroms){
    # Find centromeres
    tmp <- transScore %>% filter(chrom == chr) %>% select(control) %>% is.na() %>% as.vector() %>% rle()
    Ci <- match(max(tmp$lengths[tmp$values == T]), tmp$lengths)
    bins <- transScore[transScore$chrom == chr, 'bin']
    start <- bins[max(cumsum(tmp$lengths)[Ci-1], 0) + 1] # for those cases where Ci == 1, it will set the start as the first bin
    end <- bins[cumsum(tmp$lengths)[Ci]]
    
    # Get distances to centromeres
    transScore[with(transScore, chrom == chr & bin >= start & bin <= end), 'cenBin'] <- 0
    if(combine.arms){
      transScore[with(transScore, chrom == chr & bin < start), 'cenBin'] <- start - bins[bins < start]
    }else{
      transScore[with(transScore, chrom == chr & bin < start), 'cenBin'] <- bins[bins < start] - start
    }
    transScore[with(transScore, chrom == chr & bin > end), 'cenBin'] <- bins[bins > end] - end
  }
  
  if(plot.type == 'pileup'){
    plots_pileup <- list()
    controlIdx <- match(control, colnames(transScore))
    samples <- which(!colnames(transScore) %in% c('chrom', 'bin', 'cenBin', control))
    for(sample in samples){
      sample <- transScore[, sample]
      trend <- aggregate(sample ~ cenBin, transScore, median)
      trend$SD <- aggregate(sample ~ cenBin, transScore, sd)[,2]
      
      plot_logfc <- 
        ggplot(data = gather(transScore, sample, trans.ratio, c(samples, controlIdx) , factor_key = T)) +
        geom_point(aes(x = cenBin, y = trans.ratio, colour = sample, fill = sample), size = 0.3) + #TODO:make it reusable (names)
        #ylim(-0.5, 0.5) +
        scale_colour_manual(name = 'Cell line', values = c('grey', 'black')) +
        scale_fill_manual(name = 'Cell line', values = c('grey', 'black')) +
        geom_vline(xintercept = 0, color='black') +
        scale_x_continuous(labels = function(x) format(x/1e6)) +
        labs(x = 'Distance from centromere (Mb)', y = sprintf('Log2 fold-change of\n the proportion of %s contacts', mode)) +
        geom_line(data = trend, aes(x = cenBin, y = sample), color = 'red', size = 0.8) +
        #geom_smooth(data = gather(transScore, sample, trans.ratio, samples , factor_key = T),
        #  aes(x = cenBin, y = trans.ratio)) +
        theme_bw() +
        theme(aspect.ratio = 0.8)
      
      plots_pileup <- c(plots_pileup, list(plot_logfc))
    }
    
    ggarrange(plotlist = plots_pileup, ncol = length(samples), nrow = 1, labels = colnames(transScore)[samples])
    
  }else if (plot.type == 'perchrom'){
    transScore2 <- data.frame()
    plots_perchrom <- list()
    samples <- which(!colnames(transScore) %in% c('chrom', 'bin', 'cenBin', control))
    for(sample in samples){
      for(chr in chroms){
        dat <- transScore[transScore$chrom == chr, ]
        tmp <- aggregate(dat[, sample] ~ cenBin, dat, mean)
        names(tmp)[2] <- 'log_fc'
        tmp$chrom <- chr
        transScore2 <- rbind(transScore2, tmp)
      }
      
      plt <- ggplot(transScore2 %>% mutate(chrom = factor(chrom, levels = chroms), 
                                                            chrom = factor(chrom, levels = rev(levels(chrom)))), 
                                     aes(x = cenBin, y = chrom)) +
        geom_tile(aes(fill = log_fc, colour = log_fc)) +
        scale_fill_gradient2(low='blue', mid='white', high='red', na.value = 'grey', midpoint = 0) +
        scale_colour_gradient2(low='blue', mid='white', high='red', na.value = 'grey', midpoint = 0) +
        geom_vline(xintercept = 0, color='black') +
        labs(x = 'Distance from centromere (Mb)', y = 'Chromosome', fill = paste0('log2(', colnames(transScore)[sample], '/', control, ')')) +
        scale_x_continuous(labels = function(x) format(x/1e6)) +
        scale_y_discrete(labels = function(x) format(sub('chr', '', x))) +
        guides(colour=F) +
        GENOVA:::GENOVA_THEME() +
        theme(aspect.ratio = 0.8)
      plots_perchrom <- c(plots_perchrom, list(plt))
    }

    ggarrange(plotlist = plots_perchrom, ncol = length(samples), nrow=1, labels = colnames(transScore)[samples])

  }else{
    message('Please select plot.type either pileup or perchrom')
  }
}
