def example1():
    example_1 = '''
    n <- 50
    N <- 100
    # N <- 2*n^(1/4)*log(n)
    X <- seq(1/N,1,1/N)
    X.band <- seq(1/(N-1), 1, 1/(N-1))
    Y <- matrix(0,n,N)
    mu.true <- 5*X + 4*sin(2*pi*(X-0.5))
    psi1 <- -sqrt(2)*cos(2*pi*(X-0.5))
    psi2 <- sqrt(2) * sin(4*pi*(X-0.5))
    sig <- 0.1
    lam1 <- 2
    lam2 <- 1
    xi1 <- rnorm(n, 0, sqrt(lam1))
    xi2 <- rnorm(n, 0, sqrt(lam2))
    for(i in 1:n){
      # Y[i,]=10*sin(X)+(rnorm(N, 0, sig))
      Y[i,] <- mu.true + xi1[i]*psi1 + xi2[i]*psi2 + (rnorm(N, 0, sig))
    }
    '''
    return example_1

def example2():
    example_2 = '''
    n <- 100
    N <- 2*n^(1/4)*log(n)
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 10 + sin(2*pi*(X-1/2))
    psi1 <- -sqrt(2)*cos(pi*(X-1/2))
    psi2 <- sqrt(2)*sin(pi*(X-1/2))
    lam1 <- 2
    lam2 <- 0.5
    sig <- 0.3  # sig <- 0.5
    xi1 <- rnorm(n, 0, sqrt(lam1))
    xi2 <- rnorm(n, 0, sqrt(lam2))
    Y <- matrix(0,n,N)
    for(i in 1:n){
      # Y[i,]=10*sin(X)+(rnorm(N, 0, sig))
      Y[i,] <- mu1.true + xi1[i]*psi1 + xi2[i]*psi2 + (rnorm(N, 0, sig))
    }
    '''
    return example_2

def example3():
    example_3 = '''
    na <- 80
    nb <- 160
    N <- 100
    X <- seq(1/N,1,1/N)
    X.band <- seq(0.01, 1, 0.01)
    mu1.true <- 10 + sin(2*pi*(X-1/2))
    diff.mu <- 0.6 * X
    # diff.mu <- 0.7 * sin(X)
    # diff.mu <- 0
    mu2.true <- mu1.true + diff.mu
    psi1 <- -sqrt(2)*cos(pi*(X-1/2))
    psi2 <- sqrt(2)*sin(pi*(X-1/2))
    lam1 <- 2
    lam2 <- 0.5
    sig <- 0.3  # sig <- 0.5
    xi1.all <- list(rnorm(na,mean=0,sd=sqrt(lam1)),
                rnorm(nb,mean=0,sd=sqrt(lam1)))
    xi2.all <- list(rnorm(na,mean=0,sd=sqrt(lam2)),
                rnorm(nb,mean=0,sd=sqrt(lam2)))
    Ya <- matrix(0,na,N)
    Yb <- matrix(0,nb,N)
    for(i in 1:na){
      xi1 <- xi1.all[[1]]
      xi2 <- xi2.all[[1]]
      Ya[i,] <- mu1.true + xi1[i]*psi1 + xi2[i]*psi2 + (rnorm(N, 0, sig))
    }
    for(i in 1:nb){
      xi1 <- xi1.all[[2]]
      xi2 <- xi2.all[[2]]
      Yb[i,] <- mu2.true + xi1[i]*psi1 + xi2[i]*psi2 + (rnorm(N, 0, sig))
    }
    '''
    return example_3