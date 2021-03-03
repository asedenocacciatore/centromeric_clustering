fields.interp.surface <- function(obj, loc) {
  
  # obj is a surface or image  object like the list for contour, persp or image.
  # loc a matrix of 2 d locations -- new points to evaluate the surface.
  x <- obj$x
  y <- obj$y
  z <- obj$z
  nx <- length(x)
  ny <- length(y)
  # this clever idea for finding the intermediate coordinates at the new points
  # is from J-O Irisson
  lx <- approx(x, 1:nx, loc[, 1])$y
  ly <- approx(y, 1:ny, loc[, 2])$y
  lx1 <- floor(lx)
  ly1 <- floor(ly)
  # x and y distances between each new point and the closest grid point in the lower left hand corner.
  ex <- lx - lx1
  ey <- ly - ly1
  # fix up weights to handle the case when loc are equal to
  # last grid point.  These have been set to NA above.
  ex[lx1 == nx] <- 1
  ey[ly1 == ny] <- 1
  lx1[lx1 == nx] <- nx - 1
  ly1[ly1 == ny] <- ny - 1
  # bilinear interpolation finds simple weights based on the
  # the four corners of the grid box containing the new
  # points.
  return(z[cbind(lx1, ly1)] * (1 - ex) * 
           (1 - ey) + z[cbind(lx1 + 1, ly1)] * ex * 
           (1 - ey) + z[cbind(lx1, ly1 + 1)] * (1 - ex) * ey + 
           z[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
}

rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

resize.mat <- function(mat, ndim=dim(mat)){
  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))
  # interpolation
  ans[ncord] <- fields.interp.surface(obj, loc)
  ans[is.na(ans)] <- 0
  ans
}

expand.mat <- function(mat, ndim=dim(mat), fill=NA){
  # input object
  odim <- dim(mat)
  if(any(ndim < odim)){
    stop("New dimensions are smaller than original matrix, try resize.mat instead")
  }
  # output object
  ans <- rbind(matrix(fill, ndim[1]-odim[1], odim[2]), mat)
  ans <- cbind(matrix(fill, ndim[1], ndim[2]-odim[2]), ans)
  ans
}

select.arm <- function( mat, cp, quadrant = NULL){
  if(is.null(quadrant)){ stop("No quadrant selected") }
  
  if(quadrant==1){
    sel.i <- which(mat$x <= cp[1,2]); sel.j <- which(mat$y <= cp[1,2]);
    
  }else if(quadrant==3){	
    sel.i <- which(mat$x >= cp[1,3]); sel.j <- which(mat$y <= cp[1,2]);
    
  }else if(quadrant==2){	
    sel.i <- which(mat$x <= cp[1,2]); sel.j <- which(mat$y >= cp[1,3]);
    
  }else if(quadrant==4){	
    sel.i <- which(mat$x >= cp[1,3]); sel.j <- which(mat$y >= cp[1,3]);
  }
  
  mat$z[sel.i,sel.j]
}	

#select an arm combination of the chromosome, with two chromosomes
#that both have two arms, there are 4 possible arm combinations
#the resulting matrix is reoriented so that the centromere-centromere
#interaction is at (1,1)(bottomleft) and the telomere-telomere interaction is at
#(nrow,ncol)(topright)
select.arm.straight <- function( mat, cp, quadrant = NULL, length ){
  if(is.null(quadrant)){ stop( "No quadrant selected" ) }
  
  if(quadrant==1){
    #sel.i <- which(mat$x >= cp[1,2]-length & mat$x < cp[1,2]); sel.j <- which(mat$y >= cp[2,2]-length & mat$y < cp[2,2]);
    sel.i <- which(mat$x >= cp[1,2]-length & mat$x < cp[1,2]); sel.j <- which(mat$y >= cp[2,2]-length & mat$y < cp[2,2]);
    
  }else if(quadrant==3){	
    sel.i <- which(mat$x >= cp[1,3] & mat$x < cp[1,3] + length); sel.j <- which(mat$y >= cp[2,2]-length & mat$y < cp[2,2]);
    
  }else if(quadrant==2){	
    sel.i <- which(mat$x >= cp[1,2]-length & mat$x < cp[1,2]); sel.j <- which(mat$y >= cp[2,3] & mat$y < cp[2,3] + length);
    #sel.i <- rev(sel.i); #reorientation
  }else if(quadrant==4){	
    sel.i <- which(mat$x >= cp[1,3] & mat$x < cp[1,3] + length); sel.j <- which(mat$y >= cp[2,3] & mat$y < cp[2,3] + length);
  }
  
  mat$z[sel.i,sel.j]
}	

#select a matrix of interactions for between two chromosomes
selectTransData <- function (exp, chrom1, chrom2){
  bed <- exp$ABS
  data <- exp$ICE
  X <- bed[bed[,1]==chrom1,4]
  Y <- bed[bed[,1]==chrom2,4]
  #the order of the chromosomes matters for the analysis
  #make sure that X is smaller than Y, otherwise switch
  #them around
  if(X[1] > Y[1]){
    temp <- X
    X <- Y
    Y <- temp
    temp <- chrom1; chrom1 <- chrom2; chrom2 <- temp; #switch the chromosomes around as well
  }
  #create x and y vectors that contain the positions of the
  #entries in the matrix that we are creating
  x <- rep(X[1]:X[length(X)],        tail(Y, n=1) - Y[1] + 1)
  y <- rep(Y[1]:Y[length(Y)], each = tail(X, n=1) - X[1] + 1)
  data.sub <- data[base::list(x, y)]
  data.sub <- data.sub[!is.na(data.sub$V3)]
  #create an empty matrix, that has as many rows as the 'X' chromosome has
  #windows and as many columns as the 'Y' chromosome has windows
  mat <- matrix(0, ncol=tail(Y, n=1) - Y[1] + 1, nrow=tail(X, n=1) - X[1] + 1)
  mat[cbind(data.sub$V1-min(X)+1, data.sub$V2-min(Y)+1)] <- data.sub$V3
  x.pos <- bed[bed[,1]==chrom1,2]
  y.pos <- bed[bed[,1]==chrom2,2]
  #create a list that is compatible with the image function
  mat <- list(x=x.pos, y=y.pos, z=mat)
  mat		
}

#find the largest stretch of 0s, which is
#most likely the centromere
largest.stretch <- function( x ){
  temp <- cumsum(c(1,diff(x) - 1))
  temp2 <- rle(temp)
  x[which(temp == with(temp2, values[which.max(lengths)]))]
}	

#' Calculate a matrix of the interchromosomal interactions aligned around the centromere
#'
#' @param exp GENOVA data structure for your experiment of interest
#' @param chrom.vec a vector containing the chromosomes that should be considered in the pair-wise interchromosomal interactions
#' @param nrow dimensions of the matrix (number of columns will be the same)
#' @param leave.out two-column matrix or data.frame containing the specific chromosome combinations that need to left out
#' @param q.top top quantile of scores that should be left out of the analysis (because they are outliers)
#' @param length length of sequence around centromere that should be included in the analysis

centromere.analysis <- function( exp, chrom.vec, leave.out=NULL, q.top = 1e-5, length=40e6, centromeres = NULL){
  i.vec <- 1:(length(chrom.vec)-1)
  #empty matrix for holding the contact frequencies
  #create empty matrix
  matrix.dim = length/exp$RES
  sum.matrix <- matrix(0, ncol=matrix.dim*2, nrow=matrix.dim*2)
  num.matrix <- matrix(0, ncol=matrix.dim*2, nrow=matrix.dim*2)
  quadrant <- matrix(F, ncol=matrix.dim*2, nrow=matrix.dim*2)
  quadrant[1:matrix.dim,1:matrix.dim] <- T
  quadrant.fill <- list()
  quadrant.fill[[1]] <- quadrant
  quadrant.fill[[2]] <- quadrant[,(matrix.dim*2):1]
  quadrant.fill[[3]] <- quadrant[(matrix.dim*2):1,]
  quadrant.fill[[4]] <- quadrant[(matrix.dim*2):1,(matrix.dim*2):1]
  
  #loop over the chromosome combinations
  for( i in i.vec ){
    for(j in (i+1):(max(i.vec)+1) ){
      chrom1 <- chrom.vec[i]; chrom2 <- chrom.vec[j];
      
      #if the leave.out dataset is defined check if the chromosome combination is found
      if(!is.null(leave.out)){
        sum.val <- sum(leave.out[,1]==chrom1 & leave.out[,2]==chrom2) + sum(leave.out[,1]==chrom2 & leave.out[,2]==chrom1)
        if(sum.val > 0){
          next
        }
      }	
      cat(chrom1, "\t", chrom2, "\r")
      #select the interaction matrix between the trans chromosomes
      trans.mat <- selectTransData( exp, chrom1, chrom2 )
      
      #set the q.top highest values to this value
      th <- quantile(trans.mat$z, 1-q.top)
      trans.mat$z[trans.mat$z > th] <- th
      
      #if no centromeres have been provided
      if(is.null(centromeres)){
        #get the centromere positions empirically
        cent1 <- which(apply(trans.mat$z,1,sum)==0)
        cent2 <- which(apply(trans.mat$z,2,sum)==0)
        
        cent1 <- largest.stretch(cent1)	
        cent2 <- largest.stretch(cent2)	
        
        centromere.pos <- data.frame(chrom=c(chrom1,chrom2), start=c(min(cent1)*exp$RES-exp$RES, min(cent2)*exp$RES-exp$RES), end=c(max(cent1)*exp$RES, max(cent2)*exp$RES)) #substract one RES, to properly align the windows with centromere positions
      }else{ #use the supplied positions
        centromere.pos <- data.frame(chrom=c(chrom1,chrom2), 
                                     start=c(centromeres[centromeres[,1]==chrom1,2], centromeres[centromeres[,1]==chrom2,2] ), 
                                     end=c(centromeres[centromeres[,1]==chrom1,3], centromeres[centromeres[,1]==chrom2,3] ) )
      }
      
      #Q1,2,3,4 are the different quadrants of the chromosome-chromosome combinations
      #T------------T
      #| 1   |   2  |
      #------C------T
      #| 3   |   4  |
      #T------------T
      for( quadrant in 1:4){
        m.sub <- select.arm.straight(trans.mat, centromere.pos, quadrant, length)
        if(is.null(dim(m.sub))){ next } #acrocentric chromosomes happen
        if( nrow(m.sub) != matrix.dim || ncol(m.sub) != matrix.dim){ next } #skip if the matrix does not fit
        sum.matrix[ quadrant.fill[[quadrant]] ] <- sum.matrix[ quadrant.fill[[quadrant]] ] + m.sub
        num.matrix[ quadrant.fill[[quadrant]] ] <- num.matrix[ quadrant.fill[[quadrant]] ] + 1
      }	
      #image(sum.matrix/num.matrix)
    }
  }	
  sum.matrix/num.matrix
}

#' Calculate a matrix of the interchromosomal interactions aligned around the telomeres
#'
#' @param exp GENOVA data structure for your experiment of interest
#' @param chrom.vec a vector containing the chromosomes that should be considered in the pair-wise interchromosomal interactions
#' @param leave.out two-column matrix or data.frame containing the specific chromosome combinations that need to left out
#' @param q.top top quantile of scores that should be left out of the analysis (because they are outliers)
#' @param length length of sequence around telomere that should be included in the analysis

telomere.analysis <- function( exp, chrom.vec, leave.out=NULL, q.top = 1e-5, flank=40e6){
  i.vec <- 1:(length(chrom.vec)-1)
  #empty matrix for holding the contact frequencies
  #create empty matrix
  matrix.dim = flank/exp$RES
  sum.matrix <- matrix(0, ncol=matrix.dim+1, nrow=matrix.dim+1)
  num.matrix <- matrix(0, ncol=matrix.dim+1, nrow=matrix.dim+1)
  quadrant <- matrix(F, ncol=matrix.dim+1, nrow=matrix.dim+1)
  quadrant[1:(matrix.dim+1),1:(matrix.dim+1)] <- T
  quadrant.fill <- list()
  quadrant.fill[[1]] <- quadrant
  quadrant.fill[[2]] <- quadrant[,(matrix.dim+1):1]
  quadrant.fill[[3]] <- quadrant[,(matrix.dim+1):1]
  quadrant.fill[[4]] <- quadrant[(matrix.dim+1):1,(matrix.dim+1):1]
  
  #loop over the chromosome combinations
  for( i in i.vec ){
    for(j in (i+1):(max(i.vec)+1) ){
      chrom1 <- chrom.vec[i]; chrom2 <- chrom.vec[j];
      
      #if the leave.out dataset is defined check if the chromosome combination is found
      if(!is.null(leave.out)){
        sum.val <- sum(leave.out[,1]==chrom1 & leave.out[,2]==chrom2) + sum(leave.out[,1]==chrom2 & leave.out[,2]==chrom1)
        if(sum.val > 0){
          next
        }
      }	
      cat(chrom1, "\t", chrom2, "\r")
      #select the interaction matrix between the trans chromosomes
      trans.mat <- selectTransData( exp, chrom1, chrom2 )
      
      #set the q.top highest values to this value
      th <- quantile(trans.mat$z, 1-q.top)
      trans.mat$z[trans.mat$z > th] <- th
      
      #get the teloomere positions
      #tel1.p <- trans,mat
      #cent2 <- which(apply(trans.mat$z,2,sum)==0)
      
      
      #centromere.pos <- data.frame(chrom=c(chrom1,chrom2), start=c(min(cent1)*exp$RES-exp$RES, min(cent2)*exp$RES-exp$RES), end=c(max(cent1)*exp$RES, max(cent2)*exp$RES)) #substract one RES, to properly align the windows with centromere positions
      
      #Q1,2,3,4 are the different quadrants of the chromosome-chromosome combinations
      #T------------T
      #| 2   |   4  |
      #------C------T
      #| 1   |   3  |
      #T------------T
      
      m.sub.list <- list()
      #Q1
      m.sub <- trans.mat$z[trans.mat$x <= flank, trans.mat$y <= flank]
      sum.matrix <- sum.matrix + m.sub
      num.matrix <- num.matrix + 1
      #Q2
      m.sub <- trans.mat$z[trans.mat$x <= flank, trans.mat$y >= max(trans.mat$y) - flank]
      sum.matrix <- sum.matrix + m.sub[,(matrix.dim+1):1]
      num.matrix <- num.matrix + 1
      #Q3
      m.sub <- trans.mat$z[trans.mat$x >= max(trans.mat$x) - flank, trans.mat$y <= flank]
      sum.matrix <- sum.matrix + m.sub[(matrix.dim+1):1,]
      num.matrix <- num.matrix + 1
      #Q4
      m.sub <- trans.mat$z[trans.mat$x >= max(trans.mat$x) - flank, trans.mat$y >= max(trans.mat$y) - flank]
      sum.matrix <- sum.matrix + m.sub[(matrix.dim+1):1,(matrix.dim+1):1]
      num.matrix <- num.matrix + 1
      
      #for( quadrant in 4){
      #if( nrow(m.sub) != matrix.dim || ncol(m.sub) != matrix.dim){ next } #skip if the matrix does not fit
      #image(m.sub.list[[quadrant]], main=dim(m.sub.list[[quadrant]])); readline()
      #sum.matrix[ quadrant.fill[[quadrant]] ] <- sum.matrix[ quadrant.fill[[quadrant]] ] + m.sub.list[[quadrant]]
      #num.matrix[ quadrant.fill[[quadrant]] ] <- num.matrix[ quadrant.fill[[quadrant]] ] + 1
      #	sum.matrix <- sum.matrix + m.sub.list[[quadrant]]
      #	num.matrix <- num.matrix + 1 
      #}	
      #image(sum.matrix/num.matrix)
    }
  }	
  sum.matrix/num.matrix
  #sum.matrix
}

set.threshold <- function(m, q.th=c(0.01, 0.99)){
  zlim <- quantile(m, q.th, na.rm = T)
  m[m < zlim[[1]]] <- zlim[[1]]
  m[m > zlim[[2]]] <- zlim[[2]]
  return(m)
}

plot.pileup <- function(m, relative='centromere', diff=T, zlim=c(0, 500), signal=''){
  require(reshape2)
  require(ggplot2)
  
  if(is.null(zlim)){
    if(diff){
      m <- set.threshold(m, q.th=c(0.01, 0.99))
      z <- max(max(m), -min(m))
      zlim <- c(-z, z)
    }else{
      m <- set.threshold(m, q.th=c(0.01, 0.99))
      zlim <- c(min(m), max(m))
    }
  }else{
    m[m < zlim[[1]]] <- zlim[[1]]
    m[m > zlim[[2]]] <- zlim[[2]]
  }
  
  m.gg <- melt(m)
  colnames(m.gg) <- c('x', 'y', 'z')
  
  plt <- ggplot(m.gg, aes(x=x, y=y, fill=z)) +
    geom_tile() +
    labs(x=paste('Relative distance from', relative), y=paste('Relative distance from', relative)) +
    scale_x_discrete(limits=c(0)) +
    scale_y_discrete(limits=c(0)) +
    coord_fixed() +
    theme_bw()
  
  if(diff){
    plt + scale_fill_gradient2(low='blue', mid='white', high='red', limits = zlim, name=signal, na.value = 'white', midpoint=0)
  }else{
    higlassCol <- c('white', '#f5a623', '#d0021b', 'black')
    plt + scale_fill_gradientn(colours = higlassCol, limits = zlim, name = 'Contacts', na.value='white')
  }
}

plot.pileup.list <- function(m.list, relative='centromere', zlim=c(0, 500), diff=F, res=1e5){
  require(reshape2)
  require(ggplot2)
  
  m.gg <- data.frame()
  for (i in 1:length(m.list)){
    m <- m.list[[i]]
    if(diff){
      m <- set.threshold(m, q.th=c(0.01, 0.99))
    }
    melted <- melt(m)
    melted$sample <- names(m.list)[i]
    m.gg <- rbind(m.gg, melted)
  }
  colnames(m.gg) <- c('x', 'y', 'z', 'sample')
  m.gg$sample <- factor(m.gg$sample, levels=names(m.list))
  
  if(diff){
    z <- max(max(m.gg$z), -min(m.gg$z))
    zlim <- c(-z, z)
  }else{
    m.gg[m.gg$z < zlim[[1]], 'z'] <- zlim[[1]]
    m.gg[m.gg$z > zlim[[2]], 'z'] <- zlim[[2]]
  }
  
  plt <- ggplot(m.gg, aes(x=x, y=y, fill=z)) +
    geom_tile() +
    labs(x=paste('Relative distance from', relative, '(', res/1e3, 'Kb)'), 
         y=paste('Relative distance from', relative, '(', res/1e3, 'Kb)')) +
    coord_fixed() +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    facet_grid(.~sample) +
    GENOVA:::GENOVA_THEME()
  
  if(diff){
    plt <- plt + scale_fill_gradient2(low='blue', mid='white', high='red', limits = zlim, name='log2', na.value = 'white', midpoint=0)
  }else{
    higlassCol <- c('white', '#f5a623', '#d0021b', 'black')
    plt <- plt + scale_fill_gradientn(colours = higlassCol, limits = zlim, name = 'Contacts', na.value='white')
  }
  
  return(plt)
}

#' Calculate a matrix of intrachromosomal interactions aligned around its centromere
#' 
#' @param exp GENOVA data structure for the experiment of interest
#' @param chromsToUse a vector containing the chromosomes that should be considered
#' @param nrow dimensions of the matrix (number of columns will be the same)
#' @param q.top top quantile of scores that should be left out of the analysis (because they are outliers)
#' @param centromeres data frame containing centromere locations
#' @param bins number of proportional bins in which to divide the chromosome arm, if NULL absolute distance is used
#'
#' @return matrix with intrachromosomal interactions piled-up 
#' 
arm.analysis <- function(exp, chromsToUse, cis.only = TRUE, q.top = 1e-5, centromeres = NULL, bins = NULL){
  # Determine number of bins to consider
  if(is.null(bins)){
    # If no bins have been specified
    if(is.null(centromeres)){
      stop("Please either specify number of bins or provide centromere info")
    }else{
      # Get maximum number of bins possible
      chrom.size <- tapply(exp$ABS$V3, exp$ABS$V1, max)
      arm.size <- cbind(centromeres[1:2], 'end' = chrom.size[names(chrom.size) %in% centromeres[,1]]-centromeres[,3])
      bins <- ceiling(max(c(arm.size[,2], arm.size[,3]))/exp$RES)
    }
    mode='ABS'
  }else{
    mode='REL'
  }
  x <- c(seq(-bins, -0.5, by=1), seq(0.5, bins, by=1))
  
  # Create empty matrix for holding the contact frequencies
  sum.matrix <- matrix(0, ncol=bins*2, nrow=bins*2)
  num.matrix <- matrix(0, ncol=bins*2, nrow=bins*2)
  quadrants <- matrix(F, ncol=bins*2, nrow=bins*2)
  quadrants[1:bins,1:bins] <- T
  quadrant.fill <- list()
  quadrant.fill[[1]] <- quadrants
  quadrant.fill[[2]] <- quadrants[,(bins*2):1]
  quadrant.fill[[3]] <- quadrants[(bins*2):1,]
  quadrant.fill[[4]] <- quadrants[(bins*2):1,(bins*2):1]
  
  # Loop over the chromosomes
  for(chr in chromsToUse){
    # Select the interaction matrix for that chromosome
    mat <- select.subset(exp, chr, 0, 1e9)
    
    # Set the q.top highest to this value
    th <- quantile(mat$z, 1-q.top)
    mat$z[mat$z > th] <- th
    
    # If no centromeres were provided
    if(is.null(centromeres)){
      # Get the centromeres empirically
      cen <- which(apply(mat$z, 1, sum)==0)
      cen <- largest.stretch(cen)
      cen <- data.frame(chrom=chr, start=min(cen), end=max(cen))
    }else{
      # Get info from provided centromeres
      cen <- data.frame(chrom=chr, 
                        start=centromeres[centromeres[, 1] == chr, 2], 
                        end=centromeres[centromeres[, 1] == chr, 3])
    }
    
    # Quadrants of arm interactions
    # T------------T  T------------T
    # |  P  |  PQ  |  |  1  |   2  |
    # ------C------T  T-----C------T
    # | PQ  |   Q  |  |  3  |   4  |
    # T------------T  T------------T
    
    # For each quadrant
    for(quadrant in 1:4){
      # Get contact info (reoriented)
      m.sub <- select.arm(mat, cen, quadrant)
      if(is.null(dim(m.sub))){ next } # For acrocentric chromosomes
      #if(nrow(m.sub) != matrix.dim || ncol(m.sub) != matrix.dim){ next } #skip if the matrix does not fit
      if(mode=='REL'){
        # Scale the matrix to number of desired bins
        m.sub <- resize.mat(m.sub, ndim=c(bins, bins))
        # Add info to matrix
        sum.matrix[ quadrant.fill[[quadrant]] ] <- sum.matrix[ quadrant.fill[[quadrant]] ] + m.sub
        num.matrix[ quadrant.fill[[quadrant]] ] <- num.matrix[ quadrant.fill[[quadrant]] ] + 1
      }else if(mode=='ABS'){
        # Expand matrix to necessary size
        m.sub <- expand.mat(m.sub, ndim=c(bins, bins), 0)
        # Add info to matrix
        sum.matrix[ quadrant.fill[[quadrant]] ] <- sum.matrix[ quadrant.fill[[quadrant]] ] + m.sub
        num.matrix[ quadrant.fill[[quadrant]] ] <- num.matrix[ quadrant.fill[[quadrant]] ] + as.numeric(m.sub != 0)
      }
    }
  }
  
  res <- sum.matrix/num.matrix
  colnames(res) <- rev(x); rownames(res) <- x
  return(res)
}
