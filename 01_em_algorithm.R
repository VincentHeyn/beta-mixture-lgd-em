# 01_em_algorithm.R
library(dplyr)

# Method of Moments initialization
estBetaParams <- function(mu, var) {
  var <- max(var, 1e-6)
  mu <- pmin(pmax(mu, 1e-4), 1 - 1e-4)
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta <- alpha * (1 / mu - 1)
  return(c(alpha = max(alpha, 1e-2), beta = max(beta, 1e-2)))
}

# Mixture Cumulative Distribution Function (CDF)
pbetamix <- function(x, alpha, beta, omega) {
  sapply(x, function(val) sum(omega * pbeta(val, shape1 = alpha, shape2 = beta)))
}

# Mixture Probability Density Function (PDF)
dbetamix <- function(x, alpha, beta, omega) {
  sapply(x, function(val) sum(omega * dbeta(val, shape1 = alpha, shape2 = beta)))
}

# EM Algorithm with Information Criteria (AIC / BIC)
EM_BetaMix <- function(X, nDist = 2, maxIt = 100, epsilon = 1e-5) {
  N <- length(X)
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
    
    loglik <- sum(log(total_density))
    if (!is.na(prev_loglik) && abs(loglik - prev_loglik) < epsilon) break
    prev_loglik <- loglik
  }
  
  # Model Selection Metrics
  p <- 3 * nDist - 1 # Number of free parameters
  aic <- 2 * p - 2 * loglik
  bic <- p * log(N) - 2 * loglik
  
  return(list(
    nDist = nDist, 
    alpha = alpha, 
    beta = beta, 
    omega = omega, 
    loglik = loglik, 
    aic = aic, 
    bic = bic
  ))
}

# Calculate Regulatory Credit Risk Measures
calc_lgd_risk_measures <- function(fit, confidence = 0.999) {
  # 1. Expected LGD (Theoretical Mean)
  exp_lgd <- sum(fit$omega * (fit$alpha / (fit$alpha + fit$beta)))
  
  # 2. LGD Value at Risk (VaR) via Numerical Inversion
  target_cdf <- function(q) pbetamix(q, fit$alpha, fit$beta, fit$omega) - confidence
  var_lgd <- uniroot(target_cdf, interval = c(0.001, 0.999))$root
  
  # 3. Downturn LGD / Expected Shortfall (ES)
  tail_integrand <- function(x) x * dbetamix(x, fit$alpha, fit$beta, fit$omega)
  tail_exp <- integrate(tail_integrand, lower = var_lgd, upper = 1)$value
  es_lgd <- tail_exp / (1 - confidence)
  
  return(data.frame(
    Metric = c("Expected LGD", paste0("LGD VaR (", confidence*100, "%)"), paste0("Downturn LGD / ES (", confidence*100, "%)")),
    Value = c(exp_lgd, var_lgd, es_lgd)
  ))
}
