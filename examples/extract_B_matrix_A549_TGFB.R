#Script to extract emission matrix (B) from A549 TGF-β data
#I'm putting this script to validate the scripts

library(tidyverse)
library(MASS)


source("hiddenCOMET/R/hiddenCOMET_functions.R")



emt_data <- read.csv("hiddenCOMET/examples/real_data_examples/A549_TGFB_hidden_trajectory.csv")
cc_data <- read.csv("hiddenCOMET/examples/real_data_examples/A549_TGFB_observed_states.csv")



times_emt <- emt_data$timepoint
times_cc <- cc_data$timepoint

times <- times_emt


P <- as.matrix(emt_data[, c("Epithelial", "Hybrid", "Mesenchymal")])
P <- t(P)  
rownames(P) <- c("Epithelial", "Hybrid", "Mesenchymal")
colnames(P) <- as.character(times)


P <- apply(P, 2, function(x) {
  x <- pmax(x, 0)
  if(sum(x) == 0) rep(1/3, 3) else x/sum(x)
})



Theta <- as.matrix(cc_data[, c("G1_phase", "S_phase", "G2M_phase")])
Theta <- t(Theta)  
rownames(Theta) <- c("G1", "S", "G2M")
colnames(Theta) <- as.character(times)

Theta <- apply(Theta, 2, function(x) {
  x <- pmax(x, 0)
  if(sum(x) == 0) rep(1/3, 3) else x/sum(x)
})




observed_states <- c("G1", "S", "G2M")
hidden_states <- c("Epithelial", "Hybrid", "Mesenchymal")



B_ls <- estimate_B_ls(Theta, P)
rownames(B_ls) <- observed_states
colnames(B_ls) <- hidden_states



B_mle <- estimate_B_mle(
  Theta = Theta,
  P = P,
  observed_states = observed_states,
  hidden_states = hidden_states,
  N = 2000,       
  starts = 25,     
  method = "BFGS"
)


#validation

Theta_pred_ls <- B_ls %*% P
Theta_pred_mle <- B_mle %*% P

residuals_ls <- colSums(abs(Theta - Theta_pred_ls))
residuals_mle <- colSums(abs(Theta - Theta_pred_mle))


residual_df_ls <- data.frame(
  Timepoint = times,
  L1_residual = round(residuals_ls, 6)
)
print(residual_df_ls)

residual_df_mle <- data.frame(
  Timepoint = times,
  L1_residual = round(residuals_mle, 6)
)
print(residual_df_mle)


comparison <- data.frame(
  Observed_State = rep(observed_states, 3),
  Hidden_State = rep(hidden_states, each = 3),
  B_LS = as.vector(B_ls),
  B_MLE = as.vector(B_mle),
  Difference = as.vector(B_mle - B_ls)
)

