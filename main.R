# main.R
library(dplyr)
library(ggplot2)
library(goftest)

source("01_em_algorithm.R")

# Load Data
if (!file.exists("data/mock_lgd_data.csv")) source("00_synthetic_data.R")
data <- read.csv("data/mock_lgd_data.csv")

# Train/Test Split (80/20)
set.seed(123)
n_train <- floor(0.8 * nrow(data))
train_idx <- sample(seq_len(nrow(data)), size = n_train)
train_data <- data$lgd_diff[train_idx]
test_data <- data$lgd_diff[-train_idx]

# Model Selection: Compare K = 1, 2, 3, 4 Components via AIC & BIC
cat("--- Model Selection (AIC/BIC Evaluation) ---\n")
candidate_k <- 1:4
results_list <- list()
selection_df <- data.frame()

for (k in candidate_k) {
  fit_k <- EM_BetaMix(train_data, nDist = k)
  results_list[[k]] <- fit_k
  selection_df <- rbind(selection_df, data.frame(
    K = k,
    LogLikelihood = fit_k$loglik,
    AIC = fit_k$aic,
    BIC = fit_k$bic
  ))
}

print(selection_df)

# Pick best model based on lowest BIC
best_k <- selection_df$K[which.min(selection_df$BIC)]
best_fit <- results_list[[best_k]]
cat("\nOptimal components chosen by BIC: K =", best_k, "\n")

# Out-of-Sample Hypothesis Testing
ks_res <- ks.test(test_data, function(x) pbetamix(x, best_fit$alpha, best_fit$beta, best_fit$omega))
ad_res <- ad.test(test_data, null = function(x) pbetamix(x, best_fit$alpha, best_fit$beta, best_fit$omega))

cat("\n--- Out-of-Sample Goodness-of-Fit (Optimal Model K =", best_k, ") ---\n")
cat("Kolmogorov-Smirnov p-value:", round(ks_res$p.value, 4), "\n")
cat("Anderson-Darling p-value:  ", round(ad_res$p.value, 4), "\n")

# Risk Measures
risk_metrics <- calc_lgd_risk_measures(best_fit, confidence = 0.999)
cat("\n--- Basel III Regulatory Credit Risk Metrics ---\n")
print(risk_metrics)

# Plot
x_seq <- seq(0.001, 0.999, length.out = 300)
fit_density <- sapply(x_seq, function(x) dbetamix(x, best_fit$alpha, best_fit$beta, best_fit$omega))
emp_density <- density(test_data, from = 0, to = 1)

plot(emp_density, 
     main = paste0("Basel IRB LGD Fit & Tail Risk (K = ", best_k, " Beta Mixture)"), 
     xlab = "LGD Residuals", ylab = "Density", xlim = c(0, 1),
     ylim = c(0, max(c(fit_density, emp_density$y)) * 1.15), lwd = 2)

lines(x_seq, fit_density, col = "forestgreen", lwd = 2.5)

# Highlight 99.9% VaR cutoff line on plot
var_val <- risk_metrics$Value[2]
abline(v = var_val, col = "red", lty = 2, lwd = 2)

legend("topright", 
       legend = c("Empirical Test Data", "Optimal Beta Mixture", "99.9% LGD VaR"), 
       col = c("black", "forestgreen", "red"), lty = c(1, 1, 2), lwd = 2, bty = "n")
