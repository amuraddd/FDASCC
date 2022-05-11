######################################################################
## Example 1: (Gu et al 2014)
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
    na=nb=n#equal sample size
    sig=para[l,2]
    p.mu=p.cov=para[l,3]
    N=para[l,4]
    
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- sin(2*pi*(X-0.5)) 
    mu2.true <- 5*(X-0.6)^2 
    kappa = 2
    
    psi1.1 = -sqrt(2)*cos(pi*(X-0.5))
    psi1.2 = sqrt(2)*sin(pi*(X-0.5))
    
    psi2.1 = 1
    psi2.2 = sqrt(2)*sin(2*pi*X)
    psi2.3 = sqrt(2)*cos(2*pi*X)
    
    lam1 = c(2/5, 1)
    lam2 = c(1/10,1)
    lam3 = 1
    xi1.all <- list(rnorm(na,mean=0,sd=sqrt(lam1[1])),
                    rnorm(nb,mean=0,sd=sqrt(lam2[1])))
    xi2.all <- list(rnorm(na,mean=0,sd=sqrt(lam1[2])),
                    rnorm(nb,mean=0,sd=sqrt(lam2[2])))
    xi3.all <- list(rnorm(nb,mean=0,sd=sqrt(lam3[1])))
    
    Ya <- matrix(0,na,N)
    Yb <- matrix(0,nb,N)
    for(i in 1:na){
      xi1 <- xi1.all[[1]]
      xi2 <- xi2.all[[1]]
      Ya[i,] <- mu1.true + xi1[i]*psi1.1 + xi2[i]*psi1.2 + (rnorm(N, 0, sig))
    }
    for(i in 1:nb){
      xi1 <- xi1.all[[2]]
      xi2 <- xi2.all[[2]]
      xi3 <- xi3.all[[1]]
      Yb[i,] <- mu2.true + xi1[i]*psi2.1 + xi2[i]*psi2.2 + xi3[i]*psi2.3 + (rnorm(N, 0, sig))
    }
    
    
    Yband.true <- -sin(2*pi*(X.band-0.5))+5*(X.band-0.6)^2 
    Yband.v.true <- -2*pi*cos(2*pi*(X.band-0.5)) + 10*(X.band-0.6)
    
    out <- scc.1D(Ya=Ya,Yb=Yb, X=X,X.band=X.band, d.est=p.mu, d.cov=p.cov)
    #plot(X, out$Yhat[2,] - out$Yhat[1,], type = 'l',
         ylim = c(min(out$scc[,1,2])-1,max(out$scc[,2,2])+1))
    #lines(X.band, Yband.true, col = 2)
    #lines(X.band, out$scc[,1,2], col='blue')
    #lines(X.band, out$scc[,2,2], col='blue')
    
    colnames(result)[1:3] = c("n", "sigma", "p.mu")
    result$cov.rate[l] = mean((Yband.v.true >=  out$scc.deriv[,1,2])&(Yband.v.true <=  out$scc.deriv[,2,2]))
    result$cov.rate.v[l] = mean((Yband.true >=  out$scc[,1,2])&(Yband.true <=  out$scc[,2,2]))
    result$mse[l] = mean((mu1.true-mu2.true-as.vector(out$Yhat[1,]) + as.vector(out$Yhat[2,]))^2)
    
  }
  return(result)
}


result.1 = test.fun(n, sig, p.mu, p.cov)
