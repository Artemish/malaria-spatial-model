library(dplyr)
library(reshape2)
library(ggplot2)
library(MASS)
library(tidyr)
library(CARBayesST)
library(splines)
library(Rcpp)

load('prediction_data.Rdata')

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





final_data_lagged4[final_data_lagged4$timeid==60,'marlaria']=NA
final_data_lagged4[(nrow(final_data_lagged4)-1):nrow(final_data_lagged4),'marlaria']=NA
final_data_lagged4[nrow(final_data_lagged4),'marlaria']=NA

final_data_lagged4_aug = final_data_lagged4

frame.results <- common.frame(f, final_data_lagged4_aug, "poisson")
N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  
Y.DA <- Y     




#### Check on the rho arguments
rho <- summary['rho.S','Mean']
fix.rho.S <- FALSE   

alpha <- c(summary['rho1.T','Mean'], summary['rho2.T','Mean'])
fix.rho.T <- FALSE   

#### CAR quantities
W.quants <- common.Wcheckformat.leroux(map_mat_north)
K <- W.quants$n
N <- N.all / K
W <- W.quants$W
W.triplet <- W.quants$W.triplet
W.n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin



## Compute the blocking structure for beta     
block.temp <- common.betablock(p)
beta.beg  <- block.temp[[1]]
beta.fin <- block.temp[[2]]
n.beta.block <- block.temp[[3]]
list.block <- as.list(rep(NA, n.beta.block*2))
for(r in 1:n.beta.block)
{
  list.block[[r]] <- beta.beg[r]:beta.fin[r]-1
  list.block[[r+n.beta.block]] <- length(list.block[[r]])
}




#############################
#### Initial parameter values
#############################
beta <- summary[1:27,'Mean']
log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
res.temp <- log.Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=N.all, mean=0, sd = res.sd)
tau2 <- summary['tau2','Mean']


#### Specify matrix quantities
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)   
fitted <- exp(as.numeric(offset.mat + regression.mat + phi.mat))


###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
burnin=burn.in.car
n.sample=Ncar
thin=thinning
n.keep <- floor((n.sample - burnin)/thin)
samples.fitted <- array(NA, c(n.keep, N.all))
if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))


#### Specify the Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.05
proposal.sd.beta <- 0.01



#############################
#### Specify spatial elements
#############################
#### Spatial determinant
if(!fix.rho.S) 
{
  Wstar <- diag(apply(W,1,sum)) - W
  Wstar.eigen <- eigen(Wstar)
  Wstar.val <- Wstar.eigen$values
  det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
}else
{}



accept.all <- rep(0,6)
accept <- accept.all

#use last month for Y
if(n.miss>0)
{
  Y.DA[which.miss==0] <- tail(Y.DA[which.miss!=0],65)   
}else
{}
Y.DA.mat <- matrix(Y.DA, nrow=K, ncol=N, byrow=FALSE)

#### Create the MCMC samples
for(j in 1:n.sample) {
  
  ####################
  ## Sample from phi
  ####################
  phi.offset <- offset.mat + regression.mat
  den.offset <- rho * W.triplet.sum + 1 - rho
  temp1 <- poissonar2carupdateRW(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, alpha[1], alpha[2], rho, Y.DA.mat, proposal.sd.phi, phi.offset, den.offset)      
  phi.temp <- temp1[[1]]
  phi <- as.numeric(phi.temp)  - mean(as.numeric(phi.temp))
  phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
  
  
  lp <- as.numeric(offset.mat + regression.mat + phi.mat)
  fitted <- exp(lp)
  
  ###################
  ## Save the results
  ###################
  if(j > burnin & (j-burnin)%%thin==0)
  {
    ele <- (j - burnin) / thin
    samples.fitted[ele, ] <- fitted
  }else
  {}
  
  print(j)
}

# fitted.values <- apply(samples.fitted, 2, mean)
# fitted.values_mat = matrix(fitted.values, nrow=K, ncol=N, byrow=FALSE)
# fitted_pred = fitted.values_mat[,56]
# 
# save(fitted_pred, file = "fitted_pred_north.RData")
