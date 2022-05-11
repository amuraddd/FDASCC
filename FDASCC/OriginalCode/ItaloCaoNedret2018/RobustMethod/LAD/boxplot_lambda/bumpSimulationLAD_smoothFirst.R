########
#Confidence Band Test - Spline band----
########


source("auxCode/gervini_functions.R")
source("auxCode/auxFunctions.R")
#Libraries ----
library(splines)
library(MASS)
library(quantreg)
library(reshape2)
#Sample Generation ----

dot.prod <- function(data1, data2) {
  #TODO Need to check length of data1 and data2 are same
  N <- length(data1)
  return( sum(data1*data2) /N)
}

bumpSimLADSmoothFirst <- function(contamination, n, p = 4, Nm = 0, c = 0.5, sigma = 0.5, alpha = 0.05, strength = 10, correction = 0, debug = FALSE) { 
  #Constants ----
  
  # contamination = Proportion of dataset of outliers
  # n             = Number of samples
  # p - 1         = Degree of spline polynomials
  # Nm            = Number of knots (if 0, use formula from Cao)
  # c             = Constant used to find Nm
  # sigma         = Variance of the measrument error
  # alpha         = Confidence band
  # strength      = strengh of the outlier - dependent on type of outlier
  
  N <- floor((n^0.25) * ((log(n))^2)) #Number of observations
  
  if(Nm == 0) {
    Nm <- floor(c*(n^(1/(2*p))*log(n))) #Number of  knots
  }
  
  #Definition of Correction Factor
  if(correction == 0) {
    correction = 1.3*1.1^(log(n/30, base = 2))
  }
  #####
  
  N.test = 100
  
  #Spline pre-calculations ----
  x <- (1:N)/N
  knots<- generateKnots(p, Nm, x)
  
  #Pre-calculating the B-Splines 
  
  x.melt <- rep(x,n)
  
  X.long <- bs(x = x.melt , knots =  knots[2:(Nm+1)], degree = (p-1), intercept=TRUE, Boundary.knots=knots[c(1,Nm+2)])
  X <- bs(x = x , knots =  knots[2:(Nm+1)], degree = (p-1), intercept=TRUE, Boundary.knots=knots[c(1,Nm+2)])
  
  X.pre.calc <-  solve(t(X) %*% X) %*% t(X)
  t.coverage.test <- (1:N.test)/N.test
  
  X.coverage.test <- bs(x = t.coverage.test , knots =  knots[2:(Nm+1)], degree = (p-1), intercept=TRUE, Boundary.knots=knots[c(1,Nm+2)])
  
  real.mean <- generate.mean(N)
  mean.coverage <- generate.mean(N.test)
  
  #Beggining of simulation ----
  
  miss <- 0
  
  n.out = floor(n*contamination)
  n.real = n - n.out
  maxSim <- 500
  
  if (debug == TRUE) {
   lambda <- matrix(0, nrow = maxSim, ncol = 2)
   l2norm <- matrix(0, nrow = maxSim, ncol = 2)
   angle <-  matrix(0, nrow = maxSim, ncol = 2)
   normEig <-  matrix(0, nrow = maxSim, ncol = 2)
   phi <- generate.phi(N)
   phi[1, ]<- phi[1, ]/sqrt(2)
   phi[2, ]<- phi[2, ]/sqrt(0.5)
   realEig <- c(2,0.5)
  }
  pb <- txtProgressBar(min = 0, max = maxSim)
  
  for (a in 1 : maxSim){
    
    sample <- rbind(generateSample1(n.real, N, sigma), generateOutlierBump(numberOfSamples =  n.out, N, start = 0.9, strength = strength, sigma))
    sample.melt <- melt(t(sample))$value
    
    #Mean function Estimation ----  
    coef.mean <- as.numeric(rq(sample.melt ~ X.long - 1)$coef)
    mean.hat<- X %*% coef.mean
    
    #Covariance Estimation ----  
    sp.princomp <- f.Sp.princomp.2(t = x, x = sample, mu = as.vector(mean.hat))
    
    #Eigenvalue calculation in R is resturning negative values, instead of ~ 0, in the middle of algorithm- This amounts to returning NA for certaing scores
    # sp.princomp$scores[is.na(sp.princomp$scores)] = 0
    # eigenvalues <- apply(sp.princomp$scores, 2, function(data) hubers(data, mu = 0)$s)^2
    
    eigenvectors <- sp.princomp$pc
    eigenvectors[is.na(eigenvectors)] = 0
    eigenvectors.coef <- X.pre.calc %*% t(eigenvectors)
    eigenvectors <- t(X %*% eigenvectors.coef)
    sample.centered <- t(t(sample) - as.vector(mean.hat))
    eigenvalues <- apply(sample.centered %*% t(eigenvectors)/N, 2, function(data) hubers(data , mu=0)$s)^2
    
    eigenvalues.sum <- sum(eigenvalues, na.rm = T)
    
    for (i in 1:length(eigenvalues)){
      #Find the smallest number of eigenvalues that will make up for 95% of the total
      if(sum(eigenvalues[1:i])/eigenvalues.sum > 0.95)
      {
        k <- i
        break
      }
    }
    
    eigenvalues <- eigenvalues[1:k]
    if(debug == TRUE) {
      nEig <- min(2,k)
      for(ll in 1:nEig){
        lambda[a,ll ] <- (eigenvalues[ll]/realEig[ll]) - 1
        upNorm <- sqrt(dot.prod(eigenvectors[ll, ] - phi[ll, ], eigenvectors[ll, ] - phi[ll, ]))
        downNorm <- sqrt(dot.prod(-eigenvectors[ll, ] - phi[ll, ], -eigenvectors[ll, ] - phi[ll, ]))
        l2norm[a, ll] <- min(upNorm, downNorm)
        normEig[a, ll] <- sqrt(dot.prod(eigenvectors[ll, ], eigenvectors[ll, ]))
        upAngle <- acos(dot.prod(eigenvectors[ll, ], phi[ll, ])/normEig[a, ll])*180/pi
        downAngle <- acos(dot.prod(-eigenvectors[ll, ], phi[ll, ])/normEig[a, ll])*180/pi
        angle[a, ll] <- min(upAngle,downAngle)
      }
    }
    setTxtProgressBar(pb, a)
    
  }
  #   cat( proc.time() - timeBefore )
  close(pb)
  if(debug == TRUE){
    return(list(lambda = lambda, l2norm = l2norm, normEig = normEig,angle = angle))
  } else {
    return(miss/maxSim)  
  }
  
  
  
}

