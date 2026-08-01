# main.R
library(dplyr)
library(ggplot2)
library(goftest)

source("01_em_algorithm.R")

# 1. Load Data
if (!file.exists("data/mock_lgd_data.csv")) {
  source("00_synthetic_data.R")
}
data <- read.csv("data/mock_lgd_data.csv")

# 2. Train/Test Split (80/20)
set.seed(123)
n_train <- floor(0.8 * nrow(data))
train_idx <- sample(seq_len(nrow(data)), size = n_train)

train_data <- data$lgd_diff[train_idx]
test_data <- data$lgd_diff[-train_idx]

# 3. Fit 2-Component Beta Mixture
cat("Fitting Beta Mixture Model via EM algorithm...\n")
fit <- EM_BetaMix(train_data, nDist = 2, maxIt = 50)

cat("\n--- Estimated Parameters ---\n")
print(data.frame(Component = 1:2, Weight = fit$omega, Alpha = fit$alpha, Beta = fit$beta))

# 4. Out-of-Sample Hypothesis Testing
ks_res <- ks.test(test_data, function(x) pbetamix(x, fit$alpha, fit$beta, fit$omega))
ad_res <- ad.test(test_data, null = function(x) pbetamix(x, fit$alpha, fit$beta, fit$omega))

cat("\n--- Out-of-Sample Goodness-of-Fit Tests ---\n")
cat("Kolmogorov-Smirnov p-value:", round(ks_res$p.value, 4), "\n")
cat("Anderson-Darling p-value:  ", round(ad_res$p.value, 4), "\n")

# 5. Density Plot
x_seq <- seq(0, 1, length.out = 200)
fit_density <- sapply(x_seq, function(x) sum(fit$omega * dbeta(x, fit$alpha, fit$beta)))

plot(density(test_data), main = "Out-of-Sample LGD Density Fit", xlab = "LGD Residuals", lwd = 2)
lines(x_seq, fit_density, col = "green", lwd = 2)
legend("topright", legend = c("Empirical Test Data", "Fitted Beta Mixture"), col = c("black", "green"), lty = 1, lwd = 2)
