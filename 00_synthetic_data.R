# 00_synthetic_data.R
set.seed(42)
n <- 2500

# Simulate a 2-component Beta Mixture for LGD errors
comp1 <- rbeta(n * 0.65, shape1 = 2, shape2 = 12) # Low error cluster
comp2 <- rbeta(n * 0.35, shape1 = 6, shape2 = 3)  # High error cluster

lgd_diff <- c(comp1, comp2)
dates <- seq(as.Date("2005-01-01"), as.Date("2022-12-31"), length.out = n)

mock_data <- data.frame(
  AUSFALL_DAT = dates,
  lgd_diff = lgd_diff
)

dir.create("data", showWarnings = FALSE)
write.csv(mock_data, "data/mock_lgd_data.csv", row.names = FALSE)
cat("Synthetic dataset created at 'data/mock_lgd_data.csv'\n")
