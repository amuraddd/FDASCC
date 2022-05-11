#Two-Sample Test - Spline band
#library(ggplot2)
library(splines)



#Data Generation

#Auxiliary Functions -----
generate.phi<- function(size){
  x = (1:size)/size
  phi <- array(0.0, c(2,size))
  
  phi[1,] <- -2 * cos(pi * (x - 0.5))
  phi[2,] <-  sin(pi * (x - 0.5))
  
  return(phi)
}

generate.mean <- function(size){
  x = (1:size)/size
  
  return(10 + sqrt((x-1)^4))
}


# generate.mean <- function(size){
#   x = (1:size)/size
#   
#   return(10 + sin( 2* pi *(x-0.5)))
# }

generateSampleTDist <- function(numberOfSamples, size, sigma, dof){
  #Obtain the mean m and the phi functions
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rt(size, df = dof)
  }
  return(sample)
}



# Generate Outliers ----
####

generateOutlierPeak <- function(numberOfSamples, size, strength, sigma) {
  if(!(length(sigma) == 1 || length(sigma) == size )) {
    warning("Sigma must be a number or a vector of length 'size'.\n Sigma will be truncated.")
    sigma = sigma[1]
  }
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  ind <- round(0.05*size, digits = 0)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rnorm(size)
  }
   
  sample[ , ind] <- sample[ , ind] + strength
  return(sample)
}

generateOutlierStep <- function(numberOfSamples, size, strength, start, sigma) {
  
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  start <- round(size*runif(numberOfSamples, min = start, max = 1), digits = 0)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rnorm(size)
    sample[i, start[i]:size] <- sample[i, start[i]:size] + strength
  }
  return(sample)
}

generateOutlierBump <- function(numberOfSamples, size, strength, start = 0.5, sigma) {
  
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  # start <- floor(size*runif(numberOfSamples, min = start, max = 1))
  # end <-  pmin(start + floor(size*0.1), size)
  
  begin <- round(size*start, digits = 0)
  end <- min(round(size*(start+0.03), digits = 0), size)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rnorm(size)
    # sample[i, start[i]:end[i]] <- sample[i, start[i]:end[i]] + strength
    sample[i, begin:end] <- sample[i, begin:end] + strength
  }
  return(sample)
}

generateOutlierCauchy <- function(numberOfSamples, size, strength=0.0, start = 0.0, sigma) {
  
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  # start <- floor(size*runif(numberOfSamples, min = start, max = 1))
  # end <-  pmin(start + floor(size*0.1), size)
  
  # begin <- floor(size*start)
  # end <- min(floor(size*(start+0.03)), size)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rcauchy(size)
    # sample[i, start[i]:end[i]] <- sample[i, start[i]:end[i]] + strength
  }
  return(sample)
}
# library(rmutil)
# generateOutlierLaplace <- function(numberOfSamples, size, strength=0.0, start = 0.0, sigma) {
#   
#   m <- generate.mean(size) 
#   phi <- generate.phi(size)
#   # start <- floor(size*runif(numberOfSamples, min = start, max = 1))
#   # end <-  pmin(start + floor(size*0.1), size)
#   
#   # begin <- floor(size*start)
#   # end <- min(floor(size*(start+0.03)), size)
#   
#   sample <- array(0.0, c(numberOfSamples, size))
#   
#   for ( i in 1:numberOfSamples){
#     xsi <- array(data = rnorm(n=2), dim = dim(phi))
#     
#     sample[i,] <- m + colSums(phi*xsi) +  sigma*rlaplace(size)
#     # sample[i, start[i]:end[i]] <- sample[i, start[i]:end[i]] + strength
#   }
#   return(sample)
# }

generateOutlierBumpG <- function(numberOfSamples, size, strength, start = 0.5, sigma) {
  
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  # start <- floor(size*runif(numberOfSamples, min = start, max = 1))
  # end <-  pmin(start + floor(size*0.1), size)
  
  begin <- floor(size*start)
  end <- min(floor(size*(start+0.03)), size)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rcauchy(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi)
    # sample[i, start[i]:end[i]] <- sample[i, start[i]:end[i]] + strength
  }
  return(sample)
}

generateOutlierShape <- function(numberOfSamples, size, strength) {
  x = (1:size)/size
  
  
  m <- 10 + sin( 10* pi *(x-0.5))
  phi <- generate.phi(size)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rnorm(size) 
  }
  return(sample)
}

generateOutlierAsymmetric <- function(numberOfSamples, size, strength, sigma) {
  
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rnorm(size) + strength
  }
  return(sample)
}

generateOutlierSymmetric <- function(numberOfSamples, size, strength, sigma) {
  
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  n.up <- floor(numberOfSamples/2) + (numberOfSamples %% 2)
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi)) 
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rnorm(size) 
  }
  
  sample <- sample + c(rep(strength, n.up), rep(-1*strength, numberOfSamples-n.up))
  
  return(sample)
}



generateOutlierMixed <- function(numberOfSamples, size, strength, sigma) {
  n1 <- floor(0.5 * numberOfSamples)
  n2 <- numberOfSamples - n1
  sample1 <- generateOutlierBump(n1, size, strength, sigma = sigma)
  # sample2 <- generateOutlierSymmetric(n2, size, strength, sigma)
  sample2 <- generateOutlierBump(n1, size, strength,start = 0.8, sigma = sigma)
  return(rbind(sample1, sample2))
}



#Two-sample mean test generation----
####
generateSample1 <- function(numberOfSamples, size, sigma){
  #Obtain the mean m and the phi functions
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  
  sample <- matrix(0.0, nrow = numberOfSamples, ncol = size)
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rnorm(size)
  }
  return(sample)
}

generateSampleG <- function(numberOfSamples, size, sigma){
  #Obtain the mean m and the phi functions
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  
  sample <- matrix(0.0, nrow = numberOfSamples, ncol = size)
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rcauchy(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi)
  }
  return(sample)
}

generateSampleCauchy <- function(numberOfSamples, size, sigma){
  #Obtain the mean m and the phi functions
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  
  sample <- matrix(0.0, nrow = numberOfSamples, ncol = size)
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + colSums(phi*xsi) +  sigma*rcauchy(size)
  }
  return(sample)
}

generateDelta <- function(size, case.function){
  x<- (1:size)/size
  if(case.function == 1){
    return (2 * x)
  } else if (case.function == 2){
    return (0.7*sin(x))
  } else if (case.function == 3) {
    return (rep(0 , size))
  }
}

generateSample2 <- function(numberOfSamples, size, sigma, case.function){
  #Obtain the mean m and the phi functions
  m <- generate.mean(size) 
  phi <- generate.phi(size)
  delta <- generateDelta(size, case.function)
  
  sample <- array(0.0, c(numberOfSamples, size))
  
  for ( i in 1:numberOfSamples){
    xsi <- array(data = rnorm(n=2), dim = dim(phi))
    
    sample[i,] <- m + delta + colSums(phi*xsi) +  sigma*rnorm(size)
  }
  return(sample)
}



#Other Aux Functions ----
####
generateKnots <- function(degree, Nm, x){
  
  return( c(0, seq(1/(Nm+1), Nm/(Nm+1), 1/(Nm+1)),1) ) 
  
  # return( seq(x[1],x[length(x)],len=Nm+2) )
}