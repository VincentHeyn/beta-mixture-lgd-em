# 01_em_algorithm.R
library(dplyr)

# Beta parameter estimation using Method of Moments
estBetaParams <- function(mu, var) {
  var <- max(var, 1e-6) # Prevent division by zero
  mu <- pmin(pmax(mu, 1e-4), 1 - 1e-4)
  
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta <- alpha * (1 / mu - 1)
  
  alpha <- max(alpha, 1e-2)
  beta <- max(beta, 1e-2)
  return(c(alpha = alpha, beta = beta))
}

# Vectorized Beta Mixture CDF
pbetamix <- function(x, alpha, beta, omega) {
  sapply(x, function(val) {
    sum(omega * pbeta(val, shape1 = alpha, shape2 = beta))
  })
}

# Clean Expectation-Maximization Algorithm for Beta Mixtures
EM_BetaMix <- function(X, nDist = 2, maxIt = 50, epsilon = 1e-4) {
  N <- length(X)
  
  # Quantile-based initialization
  quantiles <- quantile(X, probs = seq(0, 1, length.out = nDist + 1))
  alpha <- numeric(nDist)
  beta <- numeric(nDist)
  omega <- numeric(nDist)
  
  for (k in 1:nDist) {
    interval <- X[X >= quantiles[k] & X <= quantiles[k + 1]]
    if (length(interval) < 2) interval <- X
    
    params <- estBetaParams(mean(interval), var(interval))
    alpha[k] <- params["alpha"]
    beta[k] <- params["beta"]
    omega[k] <- 1 / nDist
  }
  
  prev_loglik <- -Inf
  
  for (it in 1:maxIt) {
    # E-step
    gamma <- matrix(0, nrow = N, ncol = nDist)
    for (d in 1:nDist) {
      gamma[, d] <- omega[d] * dbeta(X, alpha[d], beta[d])
    }
    
    total_density <- rowSums(gamma) + 1e-12
    gamma <- gamma / total_density
    
    # M-step
    for (d in 1:nDist) {
      Nk <- sum(gamma[, d])
      omega[d] <- Nk / N
      
      mu_hat <- sum(gamma[, d] * X) / Nk
      var_hat <- sum(gamma[, d] * (X - mu_hat)^2) / Nk
      
      params <- estBetaParams(mu_hat, var_hat)
      alpha[d] <- params["alpha"]
      beta[d] <- params["beta"]
    }
    
    # Check convergence via Log-Likelihood
    loglik <- sum(log(total_density))
    if (!is.na(prev_loglik) && abs(loglik - prev_loglik) < epsilon) break
    prev_loglik <- loglik
  }
  
  return(list(alpha = alpha, beta = beta, omega = omega, loglik = loglik))
}
