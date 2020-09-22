
# Moran coefficient for vectors
MI.vec <- function(x,W,alternative="greater"){
  # convert x to a matrix and save names (if provided)
  x <- as.matrix(x)
  if(!is.null(colnames(x))) nams <- colnames(x)

  #####
  # Input
  # Checks
  #####
  if(0 %in% apply(x,2,sd)) warning("Constant term detected in x.")
  if(!(class(W) %in% c("matrix","Matrix","data.frame"))){
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if(class(W)!="matrix") W <- as.matrix(W)
  if(anyNA(x) | anyNA(W)) stop("Missing values detected")
  if(!(alternative %in% c("greater","lower", "two.sided"))){
    stop("Invalid input: 'alternative' must be either 'greater', 'lower', or 'two.sided'")
  }

  #####
  # Additional
  # Variables
  #####
  nx <- ncol(x)
  n <- nrow(W)
  df <- n-1 # 1 because only one x tested at a time
  S0 <- t(rep(1,n)) %*% W %*% rep(1,n)
  #S1 <- .5 * sum((W+W)^2)
  S1 <- sum((W*W)+(W*t(W))) # similar to moran.test
  S2 <- sum((rowSums(W)+colSums(W))^2)

  # projection matrix M
  M <- diag(n)-rep(1,n)%*%t(rep(1,n))/n
  MWM <- M %*% W %*% M

  #####
  # Moran's I
  #####
  # see Tiefelsdorf/ Boots 1995: p. 987 and Murakami et al. 2017: p. 71
  out <- data.frame(matrix(NA,nrow=nx,ncol=6))
  colnames(out) <- c("I","EI","VarI","zI","p-value","")
  for(i in 1:nx){
    # observed
    out[i,"I"] <- (n/S0) * t(x[,i]) %*% MWM %*% x[,i] / crossprod(x[,i],M) %*% x[,i]
    # expected
    out[i,"EI"] <- -1/(n-1)
    # variance - normality assumption
    out[i,"VarI"] <- ((n^2*S1 - n*S2 + 3*S0^2) / (S0^2 * (n^2-1))) - out[i,"EI"]^2
    # test statistic
    out[i,"zI"] <- (out[i,"I"]-out[i,"EI"])/sqrt(out[i,"VarI"])
    # p-value
    out[i,"p-value"] <- pfunc(z=out[i,"zI"],alternative=alternative)
    out[i,6] <- star(p=out[i,"p-value"])
  }
  if(!is.null(colnames(x))) rownames(out) <- nams

  #####
  # Output
  #####
  return(out)
}
