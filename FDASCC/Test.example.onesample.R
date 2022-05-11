
#source("scc1g.1D.R")
#source("scc2g.1D.R")
#source("fit.cov.1D.R")
#source("fit.mean.1D.R")
#source("scc.1D.R")
######################################################################
## Example 1: (Cao et al 2016. Simulation Studies)
n <- c(200, 300, 500, 800, 1200)
#N <- 4*floor(n^0.3*log(n))
sig <- c(0.1, 0.2) 
p.mu = c(2, 4) #linear and cubic
p.cov = c(2, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- sin(2*pi*(X-1/2))
    kappa = 1000
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (1/4)^(floor((1:kappa)/2))
    
    for(k in 1:500){
      psi[,2*k-1] = sqrt(2)*(cos(2*k*pi*X))
      psi[,2*k] = sqrt(2)*(sin(2*k*pi*X))
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    
    mu1.true.v <- 2*pi*cos(2*pi*(X - 1/2))
    Yband.true <- sin(2*pi*(X.band-1/2))
    Yband.v.true <- 2*pi*cos(2*pi*(X.band - 1/2))
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    

    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result1 = test.fun(n, sig, p.mu, p.cov)


######################################################################
## Example 2: (Cao et al 2012; JSPI)
n <- c(30, 50, 100, 200)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)

  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 4*X + 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))
    kappa = 8
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (2)^(-((1:kappa)-1))
    
    for(k in 1:8){
      psi[,k] = sqrt(2)*(sin(k*pi*X))
      
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    
    mu1.true.v <- 4 + 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*(0.5-X)/(0.01)
    Yband.true <- 4*X.band + 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))
    Yband.v.true <- 4 + 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*(0.5-X.band)/(0.01)
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result2 = test.fun(n, sig, p.mu, p.cov)

######################################################################
## Example 3: (Cao 2014. Simulation Studies)
n <- c(200, 300, 500, 800, 1200)
#N <- 4*floor(n^0.3*log(n))
sig <- c(0.1, 0.2) 
p.mu = c(2, 4)
p.cov = c(2, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  para$N = para$Var1
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 5*X + 4*sin(2*pi*(X-0.5)) 
    kappa = 2
    
    psi1 = sqrt(2)*cos(pi*X)
    psi2 = sqrt(2)*sin(pi*X)
    lam1 = 2/3
    lam2 = 1/2
    xi1 <- rnorm(n, 0, sqrt(lam1))
    xi2 <- rnorm(n, 0, sqrt(lam2))
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      Y[i,] <- mu1.true + xi1[i]*psi1 + xi2[i]*psi2 + (rnorm(N, 0, sig))
    }
    
    mu1.true.v <-  5 + 8*pi*cos(2*pi*(X-0.5)) 
    Yband.true <- 5*X.band + 4*sin(2*pi*(X.band-0.5)) 
    Yband.v.true <- 5 + 8*pi*cos(2*pi*(X.band-0.5)) 
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result3 = test.fun(n, sig, p.mu, p.cov)


######################################################################
## Example 4: (Combination of "Cao et al 2012; JSPI" and "Cao 2016" )
## same eigen functions in "Cao et al 2012; JSPI"
n <- c(30, 50, 100, 200, 300, 500)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 4*X + 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2)) + sin(2*pi*(X-1/2))
    kappa = 8
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (2)^(-((1:kappa)-1))
    
    for(k in 1:8){
      psi[,k] = sqrt(2)*(sin(k*pi*X))
      
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    
    mu1.true.v <- 4 + 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*(0.5-X)/(0.01)+2*pi*cos(2*pi*(X - 1/2))
    Yband.true <- 4*X.band + 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2)) + sin(2*pi*(X.band-1/2))
    Yband.v.true <- 4 + 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*(0.5-X.band)/(0.01)+2*pi*cos(2*pi*(X.band - 1/2))
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result4 = test.fun(n, sig, p.mu, p.cov)

######################################################################
## Example 5: (Combination of "Cao et al 2012; JSPI" and "Cao 2016" )
## same eigen functions in "Cao 2016"
n <- c(30, 50, 100, 200, 300, 500)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 4*X + 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2)) + sin(2*pi*(X-1/2))
    kappa = 1000
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (1/4)^(floor((1:kappa)/2))
    
    for(k in 1:500){
      psi[,2*k-1] = sqrt(2)*(cos(2*k*pi*X))
      psi[,2*k] = sqrt(2)*(sin(2*k*pi*X))
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    mu1.true.v <- 4 + 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*(0.5-X)/(0.01)+2*pi*cos(2*pi*(X - 1/2))
    Yband.true <- 4*X.band + 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2)) + sin(2*pi*(X.band-1/2))
    Yband.v.true <- 4 + 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*(0.5-X.band)/(0.01)+2*pi*cos(2*pi*(X.band - 1/2))
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result5 = test.fun(n, sig, p.mu, p.cov)

######################################################################
## Example 6: (Combination of "Cao et al 2012; JSPI" and "Cao 2016" )
## product of exp() and sin() functions
## same eigen functions in "Cao 2012"
n <- c(30, 50, 100, 200, 300, 500)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*sin(2*pi*(X-1/2))
    kappa = 8
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (2)^(-((1:kappa)-1))
    
    for(k in 1:8){
      psi[,k] = sqrt(2)*(sin(k*pi*X))
      
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    mu1.true.v <- 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*(0.5-X)/(0.01)*sin(2*pi*(X-1/2))+2*pi*cos(2*pi*(X - 1/2))*1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))
    Yband.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*sin(2*pi*(X.band-1/2))
    Yband.v.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*(0.5-X.band)/(0.01)*sin(2*pi*(X.band-1/2))+2*pi*cos(2*pi*(X.band - 1/2))*1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result6 = test.fun(n, sig, p.mu, p.cov)




######################################################################
## Example 7: (Combination of "Cao et al 2012; JSPI" and "Cao 2016" )
## product of exp() and sin() functions
## same eigen functions in "Cao 2016"
n <- c(30, 50, 100, 200, 300, 500)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*sin(2*pi*(X-1/2))
    kappa = 1000
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (1/4)^(floor((1:kappa)/2))
    
    for(k in 1:500){
      psi[,2*k-1] = sqrt(2)*(cos(2*k*pi*X))
      psi[,2*k] = sqrt(2)*(sin(2*k*pi*X))
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    mu1.true.v <- 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*(0.5-X)/(0.01)*sin(2*pi*(X-1/2))+2*pi*cos(2*pi*(X - 1/2))*1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))
    Yband.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*sin(2*pi*(X.band-1/2))
    Yband.v.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*(0.5-X.band)/(0.01)*sin(2*pi*(X.band-1/2))+2*pi*cos(2*pi*(X.band - 1/2))*1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result7 = test.fun(n, sig, p.mu, p.cov)




######################################################################
## Example 8: (Combination of "Cao et al 2012; JSPI" and "Cao 2016" )
## product of exp() and sin() functions
## same eigen functions in "Cao 2012"
n <- c(30, 50, 100, 200, 300, 500)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*4*X
    kappa = 8
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (2)^(-((1:kappa)-1))
    
    for(k in 1:8){
      psi[,k] = sqrt(2)*(sin(k*pi*X))
      
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    mu1.true.v <- 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*(0.5-X)/(0.01)*4*X+4*1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))
    Yband.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*4*X.band
    Yband.v.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*(0.5-X.band)/(0.01)*4*X.band+4*1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result8 = test.fun(n, sig, p.mu, p.cov)




######################################################################
## Example 9: (Combination of "Cao et al 2012; JSPI" and "Cao 2016" )
## product of exp() and sin() functions
## same eigen functions in "Cao 2016"
n <- c(30, 50, 100, 200, 300, 500)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*4*X
    kappa = 1000
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (1/4)^(floor((1:kappa)/2))
    
    for(k in 1:500){
      psi[,2*k-1] = sqrt(2)*(cos(2*k*pi*X))
      psi[,2*k] = sqrt(2)*(sin(2*k*pi*X))
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    mu1.true.v <- 1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))*(0.5-X)/(0.01)*4*X+4*1/(sqrt(2*pi)*0.1)*exp(-(X-0.5)^2/(2*(0.1)^2))
    Yband.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*4*X.band
    Yband.v.true <- 1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))*(0.5-X.band)/(0.01)*4*X.band+4*1/(sqrt(2*pi)*0.1)*exp(-(X.band-0.5)^2/(2*(0.1)^2))
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result9 = test.fun(n, sig, p.mu, p.cov)



######################################################################
## Example 10: (Combination of "Cao et al 2012; JSPI" and "Cao 2016" )
## product of exp() and sin() functions
## same eigen functions in "Cao 2012"
n <- c(30, 50, 100, 200, 300, 500)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 4*X*sin(2*pi*(X-1/2))
    kappa = 8
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (2)^(-((1:kappa)-1))
    
    for(k in 1:8){
      psi[,k] = sqrt(2)*(sin(k*pi*X))
      
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    mu1.true.v <- 4*sin(2*pi*(X-1/2))+2*pi*cos(2*pi*(X - 1/2))*4*X
    Yband.true <- 4*X.band*sin(2*pi*(X.band-1/2))
    Yband.v.true <- 4*sin(2*pi*(X.band-1/2))+2*pi*cos(2*pi*(X.band - 1/2))*4*X.band
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result10 = test.fun(n, sig, p.mu, p.cov)




######################################################################
## Example 11: (Combination of "Cao et al 2012; JSPI" and "Cao 2016" )
## product of exp() and sin() functions
## same eigen functions in "Cao 2016"
n <- c(30, 50, 100, 200, 300, 500)
sig <- c(0.1, 0.2) 
p.mu = c(2, 3, 4)
p.cov = c(2, 3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  
  para$N = 4*floor((para$Var1)^0.3*log(para$Var1))
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 4*X*sin(2*pi*(X-1/2))
    kappa = 1000
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = (1/4)^(floor((1:kappa)/2))
    
    for(k in 1:500){
      psi[,2*k-1] = sqrt(2)*(cos(2*k*pi*X))
      psi[,2*k] = sqrt(2)*(sin(2*k*pi*X))
    }
    
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    mu1.true.v <- 4*sin(2*pi*(X-1/2))+2*pi*cos(2*pi*(X - 1/2))*4*X
    Yband.true <- 4*X.band*sin(2*pi*(X.band-1/2))
    Yband.v.true <- 4*sin(2*pi*(X.band-1/2))+2*pi*cos(2*pi*(X.band - 1/2))*4*X.band
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mse[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mse.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result11 = test.fun(n, sig, p.mu, p.cov)

######################################################################
## Example 12: (Cao&Wang 2018. Simulation Studies)
n <- c(200, 300, 500, 800, 1200)
#N <- 4*floor(n^0.3*log(n))
sig <- c(0.1, 0.2) 
p.mu = c(3, 4)
p.cov = c(3, 4)
######################################################################
## Test function: Return coverage rate, mse, 
test.fun = function(n,sig,p.mu,p.cov){
  para = result = expand.grid(n, sig, p.mu)
  para$N = para$Var1
  for(l in 1:(dim(para)[1])){
    n=para[l,1]
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 10 + sin(2*pi*(X-0.5)) 
    kappa = 4
    
    kappa = 4
    xi = matrix(rnorm(n*kappa, 0, 1), n, kappa)
    psi = matrix(0, N, kappa)
    lam = rep(NA, kappa)
    
    for(k in 1:2){
      psi[,2*k-1] = sqrt(2)*(cos(k*pi*X))
      psi[,2*k] = sqrt(2)*(sin(k*pi*X))
      lam[2*k-1] = lam[2*k] =(0.5)^(2*k)
    }
    
    Y <- matrix(0,n,N)
    for(i in 1:n){
      S = 0
      for(k in 1:kappa){
        S = S + xi[i,k]*lam[k]*psi[,k]
      }
      Y[i,] <- mu1.true + S +(rnorm(N, 0, sig)) 
    }
    
    mu1.true.v <-  2*pi*cos(2*pi*(X-0.5)) 
    Yband.true <- 10 + sin(2*pi*(X.band-0.5)) 
    Yband.v.true <- 2*pi*cos(2*pi*(X.band-0.5)) 
    
    out <- scc.1D(Ya=Y,X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mu[l] = mean((mu1.true-as.vector(out$Yhat))^2)
    result$mu.v[l] = mean((mu1.true.v-as.vector(out$Yhat.deriv))^2)
  }
  return(result)
}


result12 = test.fun(n, sig, p.mu, p.cov)


