#' Project vector onto probability simplex
#'
#' Projects a vector onto the probability simplex (sum = 1, all elements >= 0)
#' using the algorithm from Duchi et al. (2008).
#'
#' @param v Numeric vector to project
#' @return Numeric vector projected onto simplex
#' @examples
#' proj_simplex(c(0.3, 0.8, -0.1))
#' @references Duchi, J., Shalev-Shwartz, S., Singer, Y., & Chandra, T. (2008). 
#' Efficient projections onto the l1-ball for learning in high dimensions. ICML.
#' @export
proj_simplex <- function(v){
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u * seq_along(u) > (cssv - 1)))
  th  <- (cssv[rho] - 1) / rho
  w   <- pmax(v - th, 0)
  s <- sum(w)
  if(s == 0) w[] <- 1/length(w) else w <- w/s
  w
}

#' Project matrix columns onto probability simplex
#'
#' Projects each column of a matrix onto the probability simplex.
#'
#' @param B Matrix whose columns should be projected
#' @return Matrix with columns projected onto simplex
#' @examples
#' B <- matrix(runif(9), 3, 3)
#' proj_columns_simplex(B)
#' @export
proj_columns_simplex <- function(B){
  for(j in seq_len(ncol(B))) B[, j] <- proj_simplex(B[, j])
  B
}

#' Check matrix conditioning for numerical stability
#'
#' Checks if a generator matrix is well-conditioned for eigenvalue computation.
#'
#' @param G Generator matrix to check
#' @return Logical indicating if matrix is well-conditioned
#' @examples
#' G <- matrix(c(-1, 1, 0, 0, -1, 1, 0, 0, -1), 3, 3)
#' check_matrix_conditioning(G)
#' @export
check_matrix_conditioning <- function(G) {
  
  if(any(is.nan(G)) || any(is.infinite(G)) || any(is.na(G))) {
    return(FALSE)
  }
  
  tryCatch({
    eigenvals <- eigen(G, only.values = TRUE)$values
    
    #Check for very large eigenvalues or NaN/Inf values
    if(any(is.na(eigenvals)) || any(is.infinite(eigenvals)) || max(abs(eigenvals)) > 100) {
      return(FALSE)
    }
    return(TRUE)
  }, error = function(e) {
    
    return(FALSE)
  })
}

#' Build generator matrix with linear topology
#'
#' Constructs a generator matrix for a continuous-time Markov chain with linear topology
#' (states connected sequentially: 1 <-> 2 <-> 3 <-> ... <-> n).
#'
#' @param rates Vector of rate parameters (length = 2*(n-1))
#' @param states Character vector of state names
#' @return Generator matrix G with proper conditioning
#' @examples
#' rates <- c(0.1, 0.05, 0.2, 0.15)
#' states <- c("E", "H", "M")
#' G <- build_G(rates, states)
#' @details The generator matrix G satisfies:
#' \itemize{
#'   \item G[i,i] = -sum of off-diagonal elements in row i
#'   \item G[i,j] >= 0 for i != j
#'   \item Each row sums to zero
#' }
#' For 3 states (E, H, M), rates = c(muE, muM, lamE, lamM) where:
#' muE = E → H, muM = M → H, lamE = H → E, lamM = H → M
#' @export
build_G <- function(rates, states){
  n <- length(states)
  
  if(n != 3 || !identical(states, c("E", "H", "M"))) {
    stop("This function only works for 3 states: E, H, M")
  }
  if(length(rates) != 4) {
    stop("Need exactly 4 rates: c(muE, muM, lamE, lamM)")
  }
  
  # Direct parameterization matching the app:
  # rates = c(muE, muM, lamE, lamM)
  muE  <- rates[1]  # E → H
  muM  <- rates[2]  # M → H  
  lamE <- rates[3]  # H → E
  lamM <- rates[4]  # H → M
  
  G <- rbind(c(-muE,         muE,           0),
             c(lamE,  -(lamE+lamM),      lamM),
             c(0,            muM,        -muM))
  
  rownames(G) <- colnames(G) <- c("E", "H", "M")
  
  if(!check_matrix_conditioning(G)) {
    diag(G) <- diag(G) - 1e-6
  }
  
  G
}

#' Negative log-likelihood for CTMC model
#'
#' Computes the negative log-likelihood of observed data given fixed emission matrix B,
#' initial distribution p0, and generator matrix rate parameters.
#'
#' @param rates Vector of generator matrix rate parameters (optimized by optim)
#' @param theta_mat Observed state proportions matrix (n_obs × n_timepoints)
#' @param B Emission matrix (n_obs × n_hidden)
#' @param times Vector of timepoints
#' @param p0 Initial hidden state distribution
#' @param states Character vector of hidden state names
#' @param N_pseudo Pseudo-count for multinomial likelihood (default: 3000)
#' @return Negative log-likelihood value
#' @examples
#' # Example usage with sample data
#' rates <- c(0.1, 0.05, 0.2, 0.15)
#' theta_mat <- matrix(runif(12), 3, 4)
#' B <- matrix(runif(9), 3, 3)
#' times <- c(0, 1, 2, 3)
#' p0 <- c(0.8, 0.15, 0.05)
#' states <- c("E", "H", "M")
#' nll <- negloglik_G_p0fixed(rates, theta_mat, B, times, p0, states)
#' @export
negloglik_G_p0fixed <- function(rates, theta_mat, B, times, p0, states, N_pseudo = 3000){
  n_states <- length(states)
  #Normalize columns of Theta
  Theta <- apply(theta_mat, 2, function(x){
    x <- pmax(x, 0); s <- sum(x)
    if(s == 0) rep(1/n_states, n_states) else x/s
  })
  Y <- round(N_pseudo * Theta) # pseudo-counts
  
  G <- build_G(rates, states)
  t0 <- min(times)
  
 
  tryCatch({
    #similar to COMET, expm(t(G) * dt)
    St <- lapply(times, function(tt) {
      tryCatch({
        expm::expm(t(G) * (tt - t0))
      }, error = function(e) {
        # If expm fails, return identity matrix as fallback
        diag(nrow(G))
      })
    })
    
    nll <- 0
    for(i in seq_along(times)){
      pt <- St[[i]] %*% p0
      th <- pmax(B %*% pt, 1e-12)   #predicted observed states from hidden states
      nll <- nll - sum(Y[, i] * log(th))
    }
    nll
  }, error = function(e) {
    #Return large penalty if anything fails
    1e6
  })
}

#' Fit generator matrix via multi-start optimization
#'
#' Estimates the generator matrix G using multi-start optimization with either
#' Nelder-Mead or L-BFGS-B methods.
#'
#' @param Theta Observed state proportions matrix (n_obs × n_timepoints)
#' @param B Emission matrix (n_obs × n_hidden)
#' @param times Vector of timepoints
#' @param p0 Initial hidden state distribution
#' @param states Character vector of hidden state names
#' @param starts Number of random restarts (default: 30)
#' @param N_pseudo Pseudo-count for multinomial likelihood (default: 3000)
#' @param method Optimization method: "L-BFGS-B" (default: "L-BFGS-B")
#' @return List containing:
#' \itemize{
#'   \item G: Estimated generator matrix
#'   \item p0: Initial distribution (normalized)
#'   \item P: Hidden state trajectory over time
#'   \item par: Optimal rate parameters
#'   \item value: Final negative log-likelihood value
#' }
#' @examples
#' # Example with sample data
#' Theta <- matrix(runif(12), 3, 4)
#' B <- matrix(runif(9), 3, 3)
#' times <- c(0, 1, 2, 3)
#' p0 <- c(0.8, 0.15, 0.05)
#' states <- c("E", "H", "M")
#' result <- fit_G_given_p0(Theta, B, times, p0, states)
#' @export
fit_G_given_p0 <- function(Theta, B, times, p0, states, starts = 30, N_pseudo = 3000, method = "L-BFGS-B"){
  p0 <- proj_simplex(as.numeric(p0))
  best <- list(value = Inf, par = NULL)
  
  n <- length(states)
  n_params <- 2 * (n - 1)  # Linear topology: 2*(n-1) parameters
  
        for(s in seq_len(starts)){
          # Direct parameterization in rate space [0.01, 0.8]
          init <- runif(n_params, 0.01, 0.8)
          
          fit  <- optim(
            par = init,
            fn  = negloglik_G_p0fixed,
            method = "L-BFGS-B",
            lower = rep(1e-5, n_params),
            upper = rep(5, n_params),
            theta_mat = Theta, B = B, times = times, p0 = p0, states = states, N_pseudo = N_pseudo,
            control = list(maxit = 5000)
          )
          if(fit$value < best$value && fit$value < 1e5) best <- fit  # avoid failed optimizations
        }
  G_hat <- build_G(best$par, states)
  
  t0 <- min(times)
  St <- lapply(times, function(tt) {
    tryCatch({
      expm::expm(t(G_hat) * (tt - t0))
    }, error = function(e) {
      diag(nrow(G_hat))
    })
  })
  P_hat <- sapply(St, function(S) as.numeric(S %*% p0))
  P_hat <- apply(P_hat, 2, function(v){ 
    v[v < 0] <- 0
    s <- sum(v)
    if(s == 0 || is.na(s) || is.infinite(s)) {
      cat("Warning: Trajectory sum is", s, "- using uniform distribution\n")
      rep(1/length(v), length(v))  # If sum is 0 or invalid, use uniform distribution
    } else {
      v/s
    }
  })
  rownames(P_hat) <- states; colnames(P_hat) <- as.character(times)
  
  # Verify trajectories sum to 1
  col_sums <- colSums(P_hat)
  if(any(abs(col_sums - 1) > 1e-6)) {
    cat("Warning: Some trajectory columns don't sum to 1:", col_sums, "\n")
  }
  
  list(G = G_hat, p0 = p0, P = P_hat, par = best$par, value = best$value)
}

#' Softmax function for probability normalization
#'
#' Converts a vector of logits to probabilities using softmax.
#'
#' @param z Vector of logits
#' @return Vector of probabilities (sums to 1)
#' @examples
#' softmax_n(c(1, 2, 3))
#' @export
softmax_n <- function(z) {
  u <- exp(z - max(z))
  u / sum(u)
}

#' Decode parameter vector into emission matrix
#'
#' Converts unconstrained parameters to column-stochastic emission matrix B.
#'
#' @param par Parameter vector of length n_obs * n_hidden
#' @param observed_states Character vector of observed state names
#' @param hidden_states Character vector of hidden state names
#' @return Column-stochastic emission matrix B (n_obs × n_hidden)
#' @examples
#' par <- rnorm(9)
#' observed_states <- c("G1", "S", "G2M")
#' hidden_states <- c("E", "H", "M")
#' B <- decode_B(par, observed_states, hidden_states)
#' @export
decode_B <- function(par, observed_states, hidden_states) {
  n_obs <- length(observed_states)
  n_hidden <- length(hidden_states)
  B <- matrix(NA_real_, n_obs, n_hidden)
  for (j in 1:n_hidden) {
    z <- par[(n_obs*(j-1)+1):(n_obs*j)]
    B[, j] <- softmax_n(z)
  }
  rownames(B) <- observed_states
  colnames(B) <- hidden_states
  B
}

#' Estimate emission matrix using least squares
#'
#' Estimates the emission matrix B using least squares regression:
#' Theta = B * P, solving for B.
#'
#' @param Theta Observed state proportions matrix (n_obs × n_timepoints)
#' @param P Hidden state proportions matrix (n_hidden × n_timepoints)
#' @return Estimated emission matrix B (n_obs × n_hidden)
#' @details Solves B = Theta * P^T * (P * P^T)^(-1). Uses pseudoinverse
#' if P * P^T is singular.
#' @examples
#' Theta <- matrix(runif(12), 3, 4)
#' P <- matrix(runif(12), 3, 4)
#' B_est <- estimate_B_ls(Theta, P)
#' @export
estimate_B_ls <- function(Theta, P){

  #Solve: Theta = B * P  =>  B = Theta * P^T * (P * P^T)^(-1)
  
  #Debug output
  cat("LS: Theta dims:", dim(Theta), "P dims:", dim(P), "\n")
  
  #Check matrices are in correct orientation
  if(ncol(Theta) != ncol(P)) {
    stop("Theta and P must have the same number of timepoints (columns)")
  }
  
  #Is P * P^T is invertible?
  ppt <- P %*% t(P)
  cat("P * P^T dimensions:", dim(ppt), "\n")
  cat("P * P^T determinant:", det(ppt), "\n")
  
  #Use pseudoinverse if matrix is singular
  if(abs(det(ppt)) < 1e-10) {
    cat("Using pseudoinverse for singular matrix\n")
    B <- Theta %*% t(P) %*% MASS::ginv(ppt)
  } else {
    B <- Theta %*% t(P) %*% solve(ppt)
  }
  
  cat("B matrix dimensions after calculation:", dim(B), "\n")
  proj_columns_simplex(B)
}

#' Estimate emission matrix using maximum likelihood
#'
#' Estimates the emission matrix B using maximum likelihood estimation
#' with multinomial likelihood.
#'
#' @param Theta Observed state proportions matrix (n_obs × n_timepoints)
#' @param P Hidden state proportions matrix (n_hidden × n_timepoints)
#' @param observed_states Character vector of observed state names
#' @param hidden_states Character vector of hidden state names
#' @param N Pseudo-count for multinomial likelihood (default: 2000)
#' @param starts Number of random restarts (default: 20)
#' @param method Optimization method: "BFGS" (default: "BFGS")
#' @return Estimated emission matrix B (n_obs × n_hidden)
#' @examples
#' Theta <- matrix(runif(12), 3, 4)
#' P <- matrix(runif(12), 3, 4)
#' observed_states <- c("G1", "S", "G2M")
#' hidden_states <- c("E", "H", "M")
#' B_est <- estimate_B_mle(Theta, P, observed_states, hidden_states)
#' @export
estimate_B_mle <- function(Theta, P, observed_states, hidden_states, N = 2000L, starts = 20L, method = "BFGS") {
  
  stopifnot(ncol(Theta) == ncol(P))
  n_obs <- length(observed_states)
  n_hidden <- length(hidden_states)
  Tn <- ncol(Theta)
  

  Theta <- apply(Theta, 2, function(x){ x <- pmax(x,0); s <- sum(x); if (s==0) rep(1/n_obs,n_obs) else x/s })
  P     <- apply(P,     2, function(x){ x <- pmax(x,0); s <- sum(x); if (s==0) rep(1/n_hidden,n_hidden) else x/s })
  

  if (length(N) == 1) N <- rep(N, Tn)
  Y <- sweep(Theta, 2, N, `*`) 
  
  #negative log-likelihood in terms of unconstrained logits
  nll <- function(par) {
    B <- decode_B(par, observed_states, hidden_states)  
    Th_hat <- B %*% P                                  
    -sum(Y * log(pmax(Th_hat, 1e-12)))                
  }
  
  #multi-start optimization
  n_params <- n_obs * n_hidden
  best <- list(value = Inf, par = NULL)
  for (s in seq_len(starts)) {
    par0 <- rnorm(n_params, 0, 0.2)  
    fit  <- optim(par0, nll, method = method,
                  control = list(maxit = 2000, reltol = 1e-10))
    if (fit$value < best$value) best <- fit
  }
  
  decode_B(best$par, observed_states, hidden_states)
}

#' Pretty round with trailing zeros
#'
#' Formats numbers with specified decimal places and trailing zeros.
#'
#' @param x Numeric vector to format
#' @param d Number of decimal places (default: 4)
#' @return Character vector of formatted numbers
#' @examples
#' rnd(c(0.123456, 0.1, 0.9999), d = 3)
#' @export
rnd <- function(x, d=4){ 
  format(round(x, d), nsmall = d) 
}

