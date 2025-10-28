# example_usage.R
# Example usage of hiddenCOMET package functions

# Load required libraries
library(expm)
library(MASS)

# Source the functions (or load the package once installed)
source("hiddenCOMET_functions.R")

# ============================================================================
# Example 1: Estimate Emission Matrix from Time-Course Data
# ============================================================================

# Simulate observed states data (4 states, 7 timepoints)
set.seed(123)
n_obs <- 4
n_timepoints <- 7
times <- c(0, 0.5, 1, 2, 3, 5, 7)

# Create Theta matrix (observed states over time)
Theta <- matrix(0, n_obs, n_timepoints)
Theta[1,] <- c(0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15)  # state1 decreases
Theta[2,] <- c(0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65)  # state2 increases
Theta[3,] <- rep(0.15, n_timepoints)  # state3 constant
Theta[4,] <- rep(0.05, n_timepoints)  # state4 constant

rownames(Theta) <- paste0("state", 1:n_obs)
colnames(Theta) <- as.character(times)

# Create P matrix (hidden states over time)
n_hidden <- 3
P <- matrix(0, n_hidden, n_timepoints)
P[1,] <- c(0.80, 0.70, 0.60, 0.45, 0.35, 0.25, 0.20)  # Epithelial decreases
P[2,] <- c(0.15, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)  # Hybrid increases
P[3,] <- c(0.05, 0.05, 0.10, 0.20, 0.25, 0.30, 0.30)  # Mesenchymal increases

rownames(P) <- c("Epithelial", "Hybrid", "Mesenchymal")
colnames(P) <- as.character(times)

# Estimate emission matrix using least squares
cat("=== Estimating Emission Matrix (Least Squares) ===\n")
B_ls <- estimate_B_ls(Theta, P)
print(round(B_ls, 3))

# Estimate emission matrix using maximum likelihood
cat("\n=== Estimating Emission Matrix (Maximum Likelihood) ===\n")
B_mle <- estimate_B_mle(Theta, P, 
                       observed_states = rownames(Theta),
                       hidden_states = rownames(P))
print(round(B_mle, 3))

# ============================================================================
# Example 2: Fit Generator Matrix
# ============================================================================

cat("\n=== Fitting Generator Matrix ===\n")

# Define initial distribution
p0 <- c(0.80, 0.15, 0.05)
names(p0) <- c("Epithelial", "Hybrid", "Mesenchymal")

# Fit generator matrix
result <- fit_G_given_p0(Theta = Theta, 
                        B = B_ls, 
                        times = times, 
                        p0 = p0, 
                        states = rownames(P),
                        starts = 10,  # Reduced for example
                        method = "Nelder-Mead")

# Display results
cat("Generator Matrix G:\n")
print(round(result$G, 4))

cat("\nInitial Distribution p0:\n")
print(round(result$p0, 3))

cat("\nHidden State Trajectory P(t):\n")
print(round(result$P, 3))

cat("\nTransition Rate Parameters (log-scale):\n")
print(round(result$par, 4))

cat("\nFinal Negative Log-Likelihood:\n")
print(result$value)

# ============================================================================
# Example 3: Model Validation
# ============================================================================

cat("\n=== Model Validation ===\n")

# Predict observed states using fitted model
Theta_pred <- B_ls %*% result$P

cat("Predicted vs Observed Theta (first 3 timepoints):\n")
comparison <- data.frame(
  Timepoint = times[1:3],
  Observed_state1 = Theta[1, 1:3],
  Predicted_state1 = Theta_pred[1, 1:3],
  Observed_state2 = Theta[2, 1:3],
  Predicted_state2 = Theta_pred[2, 1:3]
)
print(round(comparison, 3))

# Calculate residuals
residuals <- colSums(abs(Theta - Theta_pred))
cat("\nL1 Residuals by Timepoint:\n")
print(round(residuals, 4))

# ============================================================================
# Example 4: Extract Transition Rates
# ============================================================================

cat("\n=== Transition Rate Analysis ===\n")

# Extract transition rates from generator matrix
G <- result$G
rates <- data.frame(
  Transition = character(),
  Rate_per_day = numeric(),
  stringsAsFactors = FALSE
)

for(i in 1:nrow(G)) {
  for(j in 1:ncol(G)) {
    if(i != j && G[i,j] > 0) {
      from_state <- rownames(G)[i]
      to_state <- colnames(G)[j]
      rate <- G[i,j]
      rates <- rbind(rates, data.frame(
        Transition = paste0(from_state, "→", to_state),
        Rate_per_day = rate
      ))
    }
  }
}

print(round(rates, 4))

# Calculate half-times
leave_rates <- -diag(G)
half_times <- log(2) / leave_rates
names(half_times) <- rownames(G)

cat("\nMean Leaving Rates and Half-Times:\n")
half_time_df <- data.frame(
  State = names(half_times),
  Leaving_Rate = leave_rates,
  Half_Time_Days = half_times
)
print(round(half_time_df, 4))

cat("\n=== Example Complete ===\n")

