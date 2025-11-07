#Dose dependent emissions

emt_long <- readRDS("C:/Users/annic/Downloads/Dose-MCF10A-TGFB_45.Rds")
cc_long  <- read.csv("C:/Users/annic/Downloads/MCF10A_dose_cell_cyc_fracs.csv")


dose_map <- c("1"=0, "2"=12.5, "3"=25, "4"=50, "8"=100, "9"=200, "12"=400, "13"=800)


emt_idx <- sort(unique(emt_long$time))
cc_idx  <- sort(unique(cc_long$.TimeNum))
idx     <- emt_idx[emt_idx %in% cc_idx]

emt_sub <- emt_long[emt_long$time %in% idx, c("time","variable","value")]
cc_sub  <- cc_long[ cc_long$.TimeNum %in% idx,  c(".TimeNum",".Group","Fraction","Total")]


D <- length(idx)
X <- matrix(NA_real_, nrow=3, ncol=D, dimnames=list(c("E","H","M"), paste0("dose_", idx)))
for (j in seq_len(D)) {
  tcur <- idx[j]
  sub  <- emt_sub[emt_sub$time==tcur, ]
  X["E", j] <- sub$value[sub$variable=="Epithelial"][1]
  X["H", j] <- sub$value[sub$variable=="Hybrid"][1]
  X["M", j] <- sub$value[sub$variable=="Mesenchymal"][1]
}

X <- sweep(X, 2, colSums(X), "/")


Theta <- matrix(NA_real_, nrow=3, ncol=D, dimnames=list(c("G1","S","G2M"), paste0("dose_", idx)))
Ntot  <- numeric(D)  
for (j in seq_len(D)) {
  tcur <- idx[j]
  sub  <- cc_sub[cc_sub$.TimeNum==tcur, ]
  Theta["G1", j] <- sub$Fraction[sub$.Group=="G1"][1]
  Theta["S",  j] <- sub$Fraction[sub$.Group=="S"][1]
  Theta["G2M",j] <- sub$Fraction[sub$.Group=="G2M"][1]

  if ("Total" %in% names(sub)) Ntot[j] <- unique(sub$Total)[1]
}
Theta <- sweep(Theta, 2, colSums(Theta), "/")


dose_vals <- as.numeric(dose_map[as.character(idx)])
if (any(is.na(dose_vals))) stop("Dose mapping missing for some indices. Edit 'dose_map'.")
logd <- log(pmax(dose_vals, 1e-9)) 


if (any(is.na(Ntot))) Ntot[is.na(Ntot)] <- 1

 
K <- 3   
S <- 3   
cat_names <- c("G1","S","G2M")
emt_names <- c("E","H","M")
base_k <- "G2M"  


par_to_mats <- function(par) {
  alpha <- matrix(0, nrow=K, ncol=S, dimnames=list(cat_names, emt_names))
  beta  <- matrix(0, nrow=K, ncol=S, dimnames=list(cat_names, emt_names))
  pos <- 1
  for (s in emt_names) {
    for (k in c("G1","S")) { alpha[k, s] <- par[pos]; pos <- pos+1 }
    for (k in c("G1","S")) { beta[k,  s] <- par[pos]; pos <- pos+1 }
    
  }
  list(alpha=alpha, beta=beta)
}


compute_B_at_d <- function(alpha, beta, logd_scalar) {
  B <- matrix(NA_real_, nrow=K, ncol=S, dimnames=list(cat_names, emt_names))
  for (s in emt_names) {
    logits <- c(
      "G1"  = alpha["G1", s] + beta["G1", s]*logd_scalar,
      "S"   = alpha["S",  s] + beta["S",  s]*logd_scalar,
      "G2M" = 0 
    )
    ex <- exp(logits - max(logits))  
    B[, s] <- ex / sum(ex)
  }
  B
}


nll <- function(par) {
  mats <- par_to_mats(par)
  a <- mats$alpha; b <- mats$beta
  eps <- 1e-12
  total <- 0
  for (j in seq_len(D)) {
    Bd   <- compute_B_at_d(a, b, logd[j])
    That <- Bd %*% X[, j, drop=FALSE]
    That <- pmax(That[,1], eps) 
    ##weighted multinomial (fractions * Ntot)
    y <- Theta[, j]
    w <- Ntot[j]
    total <- total - w * sum( y * log(That) )
  }
  total
}


p0 <- rep(0, S * 4)  

fit <- optim(
  par = p0,
  fn  = nll,
  method = "BFGS",
  control = list(maxit = 1000, reltol = 1e-9)
)

cat("Converged:", fit$convergence==0, "  NLL:", fit$value, "\n")


mats <- par_to_mats(fit$par)
alpha_hat <- mats$alpha
beta_hat  <- mats$beta

##B(d) for each observed dose
B_list <- vector("list", D)
names(B_list) <- paste0("dose_", idx, "_(", dose_vals, "pM)")
for (j in seq_len(D)) {
  B_list[[j]] <- compute_B_at_d(alpha_hat, beta_hat, logd[j])
}


cat("\nDose-dependent emission matrices B(d): rows=G1,S,G2M; cols=E,H,M\n")
for (nm in names(B_list)) {
  cat("\n", nm, ":\n", sep="")
  print(round(B_list[[nm]], 4))
  cat("colSums:", round(colSums(B_list[[nm]]), 4), "\n")
}


Theta_hat <- matrix(NA_real_, nrow=3, ncol=D,
                    dimnames=list(cat_names, paste0("dose_", idx)))
for (j in seq_len(D)) {
  Bd <- B_list[[j]]
  Theta_hat[, j] <- as.numeric(Bd %*% X[, j, drop=FALSE])
}


residuals <- Theta - Theta_hat
rmse <- sqrt(mean(residuals^2))
cat("\nOverall RMSE(Theta):", round(rmse, 6), "\n")

cat("\nObserved Theta (by dose):\n");  print(round(Theta, 4))
cat("\nPredicted Theta_hat (by dose):\n"); print(round(Theta_hat, 4))

cat("\nFitted alpha (log-odds vs G2M baseline), per EMT column:\n")
print(round(alpha_hat, 4))
cat("\nFitted beta (dose slopes on log-odds), per EMT column:\n")
print(round(beta_hat, 4))





stages <- rownames(Theta)           #c("G1","S","G2M")
D <- ncol(Theta)

obs <- data.frame(
  dose = rep(dose_vals, each = length(stages)),
  stage = rep(stages, times = D),
  value = as.numeric(Theta),
  type  = "Observed"
)

pred <- data.frame(
  dose = rep(dose_vals, each = length(stages)),
  stage = rep(stages, times = D),
  value = as.numeric(Theta_hat),
  type  = "Predicted"
)

df <- rbind(obs, pred)


dose_levels <- sort(unique(dose_vals))
df$dose_f <- factor(df$dose, levels = dose_levels,
                    labels = as.character(dose_levels))


suppressPackageStartupMessages(library(ggplot2))

ggplot(df, aes(x = dose_f, y = value, group = type, color = type, linetype = type)) +
  geom_line(linewidth = 1.8) +
  geom_point(size = 2, shape=8, stroke=2) +
  facet_wrap(~ stage, nrow = 1) +
  scale_color_manual(values=c("#333446", "#7F8CAA"))+
  labs(x = "TGFβ dose (pM)", y = "Fraction", color = "", linetype = "",
       title = "MCF10A - Observed vs Predicted Cell-Cycle Fractions by Dose") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = "grey85"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )->plt.A




tidyB <- NULL
for (j in seq_along(B_list)) {
  B <- B_list[[j]]
  for (s in colnames(B)) for (k in rownames(B)) {
    tidyB <- rbind(tidyB, data.frame(
      dose  = dose_vals[j],
      EMT   = s,
      Phase = k,
      p     = as.numeric(B[k, s]),
      stringsAsFactors = FALSE
    ))
  }
}
tidyB$EMT   <- factor(tidyB$EMT,   levels = c("E","H","M"))
tidyB$Phase <- factor(tidyB$Phase, levels = c("G1","S","G2M"))


labs_pM <- formatC(sort(unique(dose_vals)), format="fg", drop0trailing=TRUE)
tidyB$dose_f <- factor(tidyB$dose, levels = sort(unique(dose_vals)),
                       labels = paste0(labs_pM, " pM"))

ggplot(tidyB, aes(x = dose_f, y = p, color = Phase, group = Phase)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2, shape=8, stroke=2) +
  facet_wrap(~ EMT, nrow = 1) +
  scale_color_manual(values = c(G1="#C5B0CD", S="#A376A2", G2M="#313647")) +
  labs(x = "TGFβ dose", y = "P(Cell cycle stage | EMT state)", color = "Phase",
       title = "Dose Trajectories of Emission Probabilities by EMT State") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "grey95", color = "grey85"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )



















