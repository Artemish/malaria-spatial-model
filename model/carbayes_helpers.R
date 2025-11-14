library(Rcpp)

sourceCpp("./CARBayesST.cpp")
#### Read in and format the frame argument
common.frame <- function(formula, data, family)
{
  #### Overall formula object
  frame <- try(suppressWarnings(model.frame(formula, data=data, na.action=na.pass)), silent=TRUE)
  if(class(frame)[1]=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)
  
  
  #### Design matrix
  ## Create the matrix
  X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
  if(class(X)[1]=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
  if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
  
  n <- nrow(X)
  p <- ncol(X)
  
  ## Check for linearly related columns
  cor.X <- suppressWarnings(cor(X))
  diag(cor.X) <- 0
  if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
  if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
  if(p>1)
  {
    if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
  }else
  {
  }
  
  ## Standardise the matrix
  X.standardised <- X
  X.sd <- apply(X, 2, sd)
  X.mean <- apply(X, 2, mean)
  X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back
  
  for(j in 1:p)
  {
    if(length(table(X[ ,j]))>2)
    {
      X.indicator[j] <- 1
      X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
    }else if(length(table(X[ ,j]))==1)
    {
      X.indicator[j] <- 2
    }else
    {
      X.indicator[j] <- 0
    }
  }
  
  
  #### Offset variable
  offset <- try(model.offset(frame), silent=TRUE)
  if(class(offset)[1]=="try-error")   stop("the offset is not numeric.", call.=FALSE)
  if(is.null(offset))  offset <- rep(0,n)
  if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
  
  
  #### Response variable
  ## Create the response
  Y <- model.response(frame)
  which.miss <- as.numeric(!is.na(Y))
  n.miss <- n - sum(which.miss)
  
  
  ## Check for errors
  if(family=="binomial")
  {
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
    int.check <- n - n.miss - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
  }else if(family=="gaussian")
  {
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)    }else
    {
      if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
      int.check <- n - n.miss - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
      if(int.check > 0) stop("the response variable has non-integer values.", call.=FALSE)
      if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
    }
  
  
  #### Return the values needed
  results <- list(n=n, p=p, X=X, X.standardised=X.standardised, X.sd=X.sd, X.mean=X.mean, X.indicator=X.indicator, 
                  offset=offset, Y=Y, which.miss=which.miss, n.miss=n.miss)
  return(results)
}

#### Check the W matrix - Leroux model
common.Wcheckformat.leroux <- function(W)
{
  #### Check W is a matrix of the correct dimension
  if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
  n <- nrow(W)
  if(ncol(W)!= n) stop("W is not a square matrix.", call.=FALSE)    
  
  
  #### Check validity of inputed W matrix
  if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
  if(min(W)<0) stop("W has negative elements.", call.=FALSE)
  if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
  if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
  
  
  #### Create the triplet form
  W.triplet <- c(NA, NA, NA)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      if(W[i,j]>0)
      {
        W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
      }else{}
    }
  }
  W.triplet <- W.triplet[-1, ]     
  n.triplet <- nrow(W.triplet) 
  W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
  n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
  
  
  #### Create the start and finish points for W updating
  W.begfin <- array(NA, c(n, 2))     
  temp <- 1
  for(i in 1:n)
  {
    W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
    temp <- temp + n.neighbours[i]
  }
  
  
  #### Return the critical quantities
  results <- list(W=W, W.triplet=W.triplet, n.triplet=n.triplet, W.triplet.sum=W.triplet.sum, n.neighbours=n.neighbours, W.begfin=W.begfin, n=n)
  return(results)   
}


#### Beta blocking
common.betablock <- function(p)
{
  ## Compute the blocking structure for beta     
  blocksize.beta <- 5 
  if(blocksize.beta >= p)
  {
    n.beta.block <- 1
    beta.beg <- 1
    beta.fin <- p
  }else
  {
    n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
    remainder <- p - n.standard * blocksize.beta
    
    if(remainder==0)
    {
      beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
      beta.fin <- seq(blocksize.beta, p, blocksize.beta)
      n.beta.block <- length(beta.beg)
    }else
    {
      beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
      beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
      n.beta.block <- length(beta.beg)
    }
  }
  
  return(list(beta.beg, beta.fin, n.beta.block))
}
