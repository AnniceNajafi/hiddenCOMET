# Simple script to extract emission matrix (B) from A549 TGF-β data
# Uses hiddenCOMET package functions

library(tidyverse)
library(MASS)

# Source hiddenCOMET functions
source("../R/hiddenCOMET_functions.R")

# Read data
emt_data <- read.csv("real_data_examples/A549_TGFB_hidden_trajectory.csv")
cc_data <- read.csv("real_data_examples/A549_TGFB_observed_states.csv")

# Prepare matrices (rows = states, columns = timepoints)
P <- t(as.matrix(emt_data[, c("Epithelial", "Hybrid", "Mesenchymal")]))
Theta <- t(as.matrix(cc_data[, c("G1_phase", "S_phase", "G2M_phase")]))

# Normalize columns
P <- apply(P, 2, function(x) { x <- pmax(x, 0); x/sum(x) })
Theta <- apply(Theta, 2, function(x) { x <- pmax(x, 0); x/sum(x) })

# State names
observed_states <- c("G1", "S", "G2M")
hidden_states <- c("Epithelial", "Hybrid", "Mesenchymal")
rownames(P) <- hidden_states
rownames(Theta) <- observed_states

# Estimate B matrix using Maximum Likelihood
B <- estimate_B_mle(
  Theta = Theta,
  P = P,
  observed_states = observed_states,
  hidden_states = hidden_states,
  N = 2000,
  starts = 25,
  method = "BFGS"
)

# Display results
cat("\nEmission Matrix B (P(observed | hidden)):\n")
print(round(B, 4))

# Save result
write.csv(B, "real_data_examples/A549_TGFB_B_matrix.csv", row.names = TRUE)
cat("\nSaved to: real_data_examples/A549_TGFB_B_matrix.csv\n")

