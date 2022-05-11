#!IMPORTANT!
## Set the folllowing variable to the path to the current file
#!IMPORTANT!
generate.phi <- function(size) {
  x = (1:size)/size
  phi <- array(0, c(2, size))

  phi[1, ] <- -2 * cos(pi * (x - 0.5))
  phi[2, ] <- sin(pi * (x - 0.5))

  return(phi)
}

generate.mean <- function(size) {
  x = (1:size)/size
  # return(rep(0, size))
  return(10 + sin(2 * pi * (x - 0.5)))
}

generateOutlierCauchy <- function(numberOfSamples, size, sigma) {

  m <- generate.mean(size)
  phi <- generate.phi(size)
  # start <- round(size*runif(numberOfSamples, min = start, max = 1)) end <- pmin(start + round(size*0.1), size)

  # begin <- round(size*start) end <- min(round(size*(start+0.03)), size)

  sample <- array(0, c(numberOfSamples, size))

  for (i in 1:numberOfSamples) {
    xsi <- array(data = rnorm(n = 2), dim = dim(phi))

    sample[i, ] <- m + colSums(phi * xsi) + sigma * rcauchy(size)
    # sample[i, start[i]:end[i]] <- sample[i, start[i]:end[i]] + strength
  }
  return(sample)
}

generateSample <- function(numberOfSamples, size, sigma) {
  # Obtain the mean m and the phi functions
  m <- generate.mean(size)
  phi <- generate.phi(size)

  sample <- matrix(0, nrow = numberOfSamples, ncol = size)

  for (i in 1:numberOfSamples) {
    xsi <- array(data = rnorm(n = 2), dim = dim(phi))

    sample[i, ] <- m + colSums(phi * xsi) + sigma * rnorm(size)
  }
  return(sample)
}



generateKnots <- function(degree, Nm, x) {
  return(c(0, seq(1/(Nm + 1), Nm/(Nm + 1), 1/(Nm + 1)), 1))
}

{
  library('RSCB')
  set.seed(42)
  cat(" Simulation of Mixture-Cauchy Errors.\nSample size, n = 200\nMeasurement points, N = 220\n contamination proportion 20%\np = 4, 1 - alpha = 95%")
  invisible(readline(prompt="Press [enter] to continue:"))

  contamination = .20
  n = 50
  N = 70
  p = 4
  sigma = 0.5

  x <- 1:N/N

  n.out = floor(contamination*n)
  n.real = n - n.out

  mean_real <- generate.mean(N)

  sample <- rbind(generateSample(n.real, N, sigma), generateOutlierCauchy(numberOfSamples = n.out, size = N, sigma = sigma))
  rscb <- mRSCB(sample, k = 2.50)
  scb <- SCB(sample)

  matplot(x, t(sample), type ="l",
          col = rgb(0.8,0.8,0.8, 0.5),
          lty = 1, xlab = "x", ylab = "y", ylim = c(5,15), main = "Mixture Normal-Cauchy")

  lines(x, mean_real, lwd = 2)
  lines(x, rscb$m_hat, lty = 4, lwd = 2, col = "black")
  lines(x, scb$m_hat, lty = 4, lwd = 2, col = "red")

  polygon(c(x, rev(x)), c(rscb$upper_bound, rev(rscb$lower_bound)), col = rgb(0, 0, 0, 0.5) , border = rgb(0,0,0,0))
  polygon(c(x, rev(x)), c(scb$upper_bound, rev(scb$lower_bound)), col = rgb(1, 0, 0, 0.35) , border = rgb(0,0,0,0))
}
