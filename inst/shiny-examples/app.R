
suppressPackageStartupMessages({
  library(shiny)
  library(shinyvalidate)
  library(tidyverse)
  library(DT)
  library(expm)
  library(plotly)
  library(shinycssloaders)
  library(MASS)
})

DEFAULT_OBSERVED_STATES <- c("G1", "S", "G2M")
DEFAULT_HIDDEN_STATES <- c("E", "H", "M")
DEFAULT_TOPOLOGY <- "linear"

clean_names_safe <- function(df){
  nm <- names(df); bad <- is.na(nm) | trimws(nm) == ""
  if(any(bad)) nm[bad] <- paste0("col_", seq_len(sum(bad)))
  names(df) <- make.names(nm, unique = TRUE); df
}

proj_simplex <- function(v){
  u <- sort(v, decreasing = TRUE); cssv <- cumsum(u)
  rho <- max(which(u * seq_along(u) > (cssv - 1)))
  th  <- (cssv[rho] - 1) / rho
  w   <- pmax(v - th, 0); s <- sum(w)
  if(s == 0) w[] <- 1/length(w) else w <- w/s
  w
}

proj_columns_simplex <- function(B){
  for(j in seq_len(ncol(B))) B[, j] <- proj_simplex(B[, j])
  B
}

check_matrix_conditioning <- function(G) {
  if(any(is.nan(G)) || any(is.infinite(G)) || any(is.na(G))) {
    return(FALSE)
  }
  
  tryCatch({
    eigenvals <- eigen(G, only.values = TRUE)$values
    
    if(any(is.na(eigenvals)) || any(is.infinite(eigenvals)) || max(abs(eigenvals)) > 100) {
      return(FALSE)
    }
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

build_G <- function(eta, states){
  n <- length(states)
  G <- matrix(0, n, n)
  dimnames(G) <- list(states, states)
  
  if(n < 2) stop("Linear topology requires at least 2 states")
  if(length(eta) != 2*(n-1)) stop("Linear topology requires 2*(n-1) parameters")
  
  k <- 1
  for(i in 1:(n-1)){
    G[i, i+1] <- exp(eta[k])
    G[i, i] <- G[i, i] - exp(eta[k])
    k <- k + 1
    G[i+1, i] <- exp(eta[k])
    G[i+1, i+1] <- G[i+1, i+1] - exp(eta[k])
    k <- k + 1
  }
  
  if(!check_matrix_conditioning(G)) {
    diag(G) <- diag(G) - 1e-6
  }
  
  G
}

negloglik_G_p0fixed <- function(theta_mat, B, times, eta, p0, states, N_pseudo = 3000){
  n_states <- length(states)
  Theta <- apply(theta_mat, 2, function(x){
    x <- pmax(x, 0); s <- sum(x)
    if(s == 0) rep(1/n_states, n_states) else x/s
  })
  Y <- round(N_pseudo * Theta) # pseudo-counts
  
    G <- build_G(eta, states)
  t0 <- min(times)
  
  tryCatch({
    St <- lapply(times, function(tt) {
      tryCatch({
        expm(t(G) * (tt - t0))
      }, error = function(e) {
        diag(nrow(G))
      })
    })
    
    nll <- 0
    for(i in seq_along(times)){
      pt <- St[[i]] %*% p0
      th <- pmax(B %*% pt, 1e-12)   # predicted observed states from hidden states
      nll <- nll - sum(Y[, i] * log(th))
    }
    nll
  }, error = function(e) {
    1e6
  })
}

fit_G_given_p0 <- function(Theta, B, times, p0, states, starts = 30, N_pseudo = 3000, method = "Nelder-Mead"){
  p0 <- proj_simplex(as.numeric(p0))
  best <- list(value = Inf, par = NULL)
  
  n <- length(states)
  n_params <- 2 * (n - 1)  # Linear topology: 2*(n-1) parameters
  
        for(s in seq_len(starts)){
          if(method == "BFGS") {
            init <- log(runif(n_params, 0.05, 0.2))  # very small range for BFGS
            control_list <- list(maxit = 1000, reltol = 1e-6, abstol = 1e-6)
          } else {
            init <- log(runif(n_params, 0.01, 0.3))
            control_list <- list(maxit = 3000, reltol = 1e-8)
          }
          
          fit  <- optim(
            par = init,
            fn  = negloglik_G_p0fixed,
            method = method,
            theta_mat = Theta, B = B, times = times, p0 = p0, states = states, N_pseudo = N_pseudo,
            control = control_list
          )
          if(fit$value < best$value && fit$value < 1e5) best <- fit  # avoid failed optimizations
        }
  G_hat <- build_G(best$par, states)
  
  t0 <- min(times)
  St <- lapply(times, function(tt) {
    tryCatch({
      expm(t(G_hat) * (tt - t0))
    }, error = function(e) {
      diag(nrow(G_hat))
    })
  })
  P_hat <- sapply(St, function(S) as.numeric(S %*% p0))
  P_hat <- apply(P_hat, 2, function(v){ v[v < 0] <- 0; v/sum(v) })
  rownames(P_hat) <- states; colnames(P_hat) <- as.character(times)
  
  list(G = G_hat, p0 = p0, P = P_hat, par = best$par, value = best$value)
}

rnd <- function(x, d=4){ format(round(x, d), nsmall = d) }

softmax_n <- function(z) {
  u <- exp(z - max(z))
  u / sum(u)
}

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

estimate_B_ls <- function(Theta, P){
  
  cat("LS: Theta dims:", dim(Theta), "P dims:", dim(P), "\n")
  
  if(ncol(Theta) != ncol(P)) {
    stop("Theta and P must have the same number of timepoints (columns)")
  }
  
  ppt <- P %*% t(P)
  cat("P * P^T dimensions:", dim(ppt), "\n")
  cat("P * P^T determinant:", det(ppt), "\n")
  
  if(abs(det(ppt)) < 1e-10) {
    cat("Using pseudoinverse for singular matrix\n")
    B <- Theta %*% t(P) %*% MASS::ginv(ppt)
  } else {
    B <- Theta %*% t(P) %*% solve(ppt)
  }
  
  cat("B matrix dimensions after calculation:", dim(B), "\n")
  proj_columns_simplex(B)
}

estimate_B_mle <- function(Theta, P, observed_states, hidden_states, N = 2000L, starts = 20L, method = "BFGS") {
  stopifnot(ncol(Theta) == ncol(P))
  n_obs <- length(observed_states)
  n_hidden <- length(hidden_states)
  Tn <- ncol(Theta)
  
  Theta <- apply(Theta, 2, function(x){ x <- pmax(x,0); s <- sum(x); if (s==0) rep(1/n_obs,n_obs) else x/s })
  P     <- apply(P,     2, function(x){ x <- pmax(x,0); s <- sum(x); if (s==0) rep(1/n_hidden,n_hidden) else x/s })
  
  if (length(N) == 1) N <- rep(N, Tn)
  Y <- sweep(Theta, 2, N, `*`)  # expected counts; using fractions is fine for MLE
  
  nll <- function(par) {
    B <- decode_B(par, observed_states, hidden_states)  # n_obs x n_hidden, col-stochastic
    Th_hat <- B %*% P                                   # n_obs x T
    -sum(Y * log(pmax(Th_hat, 1e-12)))                 # multinomial loglik (up to const)
  }
  
  n_params <- n_obs * n_hidden
  best <- list(value = Inf, par = NULL)
  for (s in seq_len(starts)) {
    par0 <- rnorm(n_params, 0, 0.2)  # small random logits around equal cols
    fit  <- optim(par0, nll, method = method,
                  control = list(maxit = 2000, reltol = 1e-10))
    if (fit$value < best$value) best <- fit
  }
  
  decode_B(best$par, observed_states, hidden_states)
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Import elegant fonts */
      @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap');
      
      /* Luxurious theme base */
      body {
        background: linear-gradient(135deg, #0a0a0a 0%, #1a1520 25%, #2a1f35 50%, #1a1520 75%, #0a0a0a 100%) !important;
        color: #f5f5f5 !important;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
        font-weight: 400;
        font-size: 16px;
        line-height: 1.6;
        -webkit-font-smoothing: antialiased;
        -moz-osx-font-smoothing: grayscale;
        position: relative;
        overflow-x: hidden;
      }
      
      /* Elegant subtle pattern */
      body::before {
        content: '';
        position: fixed;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background-image: 
          radial-gradient(circle at 20% 50%, rgba(158, 202, 214, 0.03) 0%, transparent 50%),
          radial-gradient(circle at 80% 80%, rgba(158, 202, 214, 0.03) 0%, transparent 50%);
        pointer-events: none;
        z-index: -1;
      }
      
      /* Main container */
      .container-fluid {
        background: transparent !important;
      }
      
      /* Luxurious logo and title styling */
      .title-panel {
        background: linear-gradient(135deg, #1a1520, #2a1f35, #1a1520) !important;
        color: #9ECAD6 !important;
        padding: 35px 30px !important;
        border-radius: 12px !important;
        margin-bottom: 30px !important;
        box-shadow: 
          0 8px 32px rgba(0, 0, 0, 0.6),
          0 0 40px rgba(158, 202, 214, 0.15),
          inset 0 1px 0 rgba(255, 255, 255, 0.05) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        border-top: 2px solid rgba(158, 202, 214, 0.6) !important;
        position: relative !important;
        overflow: hidden !important;
      }
      
      .title-panel::before {
        content: '' !important;
        position: absolute !important;
        top: 0 !important;
        left: 0 !important;
        right: 0 !important;
        height: 1px !important;
        background: linear-gradient(90deg, transparent, rgba(158, 202, 214, 0.8), transparent) !important;
      }
      
      .logo-container {
        display: flex !important;
        align-items: center !important;
        gap: 20px !important;
      }
      
      .logo {
        font-size: 3.5rem !important;
        color: #9ECAD6 !important;
        text-shadow: 
          0 2px 4px rgba(0, 0, 0, 0.5),
          0 0 20px rgba(158, 202, 214, 0.3) !important;
        filter: drop-shadow(0 4px 8px rgba(0, 0, 0, 0.4)) !important;
      }
      
      .logo-text {
        flex: 1 !important;
      }
      
      .app-title {
        color: #9ECAD6 !important;
        font-family: 'Inter', sans-serif !important;
        font-weight: 700 !important;
        font-size: 2.8rem !important;
        margin: 0 0 10px 0 !important;
        text-shadow: 
          0 2px 4px rgba(0, 0, 0, 0.5),
          0 0 20px rgba(158, 202, 214, 0.2) !important;
        letter-spacing: 0.05em !important;
      }
      
      .app-subtitle {
        color: #9ECAD6 !important;
        font-family: 'Inter', sans-serif !important;
        font-weight: 400 !important;
        font-size: 1.25rem !important;
        margin: 0 !important;
        opacity: 0.9 !important;
        letter-spacing: 0.03em !important;
        font-style: italic !important;
      }
      
      /* Luxurious sidebar panel */
      .sidebar-panel {
        background: linear-gradient(145deg, #1a1520, #2a1f35) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        border-top: 2px solid rgba(158, 202, 214, 0.5) !important;
        border-radius: 10px !important;
        padding: 28px !important;
        box-shadow: 
          0 8px 32px rgba(0, 0, 0, 0.5),
          inset 0 1px 0 rgba(255,255,255,0.03) !important;
        backdrop-filter: blur(10px) !important;
        position: relative !important;
      }
      
      .sidebar-panel::before {
        content: '' !important;
        position: absolute !important;
        top: 0 !important;
        left: 0 !important;
        right: 0 !important;
        height: 1px !important;
        background: linear-gradient(90deg, transparent, rgba(158, 202, 214, 0.6), transparent) !important;
        border-radius: 10px 10px 0 0 !important;
      }
      
      /* Luxurious main panel */
      .main-panel {
        background: linear-gradient(145deg, #1a1520, #2a1f35) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        border-top: 2px solid rgba(158, 202, 214, 0.5) !important;
        border-radius: 10px !important;
        padding: 28px !important;
        box-shadow: 
          0 8px 32px rgba(0, 0, 0, 0.5),
          inset 0 1px 0 rgba(255,255,255,0.03) !important;
        backdrop-filter: blur(10px) !important;
        position: relative !important;
      }
      
      .main-panel::before {
        content: '' !important;
        position: absolute !important;
        top: 0 !important;
        left: 0 !important;
        right: 0 !important;
        height: 1px !important;
        background: linear-gradient(90deg, transparent, rgba(158, 202, 214, 0.6), transparent) !important;
        border-radius: 10px 10px 0 0 !important;
      }
      
      /* Luxurious form controls */
      .form-control {
        background: rgba(42, 31, 53, 0.6) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        color: #f5f5f5 !important;
        border-radius: 6px !important;
        padding: 10px 14px !important;
        font-family: 'Inter', sans-serif !important;
        font-size: 16px !important;
        font-weight: 300 !important;
        transition: all 0.3s ease !important;
        box-shadow: 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 1px 0 rgba(255,255,255,0.03) !important;
      }
      
      .form-control:focus {
        background: rgba(42, 31, 53, 0.8) !important;
        border-color: rgba(158, 202, 214, 0.6) !important;
        box-shadow: 
          0 0 0 0.2rem rgba(158, 202, 214, 0.1), 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 4px 12px rgba(158, 202, 214, 0.15) !important;
        color: #ffffff !important;
        outline: none !important;
      }
      
      .form-control:hover {
        border-color: rgba(158, 202, 214, 0.5) !important;
        box-shadow: 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 2px 8px rgba(158, 202, 214, 0.1) !important;
      }
      
      /* Luxurious file input styling */
      .file-input-wrapper {
        background-color: rgba(42, 31, 53, 0.4) !important;
        border: 1px dashed rgba(158, 202, 214, 0.4) !important;
        border-radius: 6px !important;
        padding: 15px !important;
        text-align: center !important;
        transition: all 0.3s ease !important;
      }
      
      .file-input-wrapper:hover {
        border-color: rgba(158, 202, 214, 0.6) !important;
        background-color: rgba(42, 31, 53, 0.6) !important;
      }
      
      /* Luxurious buttons */
      .btn-primary {
        background: linear-gradient(135deg, rgba(158, 202, 214, 0.9), rgba(184, 134, 11, 0.9)) !important;
        border: 1px solid rgba(158, 202, 214, 0.8) !important;
        color: #0a0a0a !important;
        font-weight: 500 !important;
        font-family: 'Inter', sans-serif !important;
        padding: 12px 28px !important;
        border-radius: 6px !important;
        box-shadow: 
          0 4px 15px rgba(158, 202, 214, 0.3), 
          inset 0 1px 0 rgba(255,255,255,0.2) !important;
        transition: all 0.3s ease !important;
        text-transform: uppercase !important;
        letter-spacing: 0.05em !important;
        font-size: 16px !important;
      }
      
      .btn-primary:hover {
        background: linear-gradient(135deg, rgba(232, 193, 75, 1), rgba(158, 202, 214, 1)) !important;
        border-color: rgba(232, 193, 75, 0.9) !important;
        transform: translateY(-2px) !important;
        box-shadow: 
          0 6px 20px rgba(158, 202, 214, 0.4), 
          inset 0 1px 0 rgba(255,255,255,0.3) !important;
        color: #000000 !important;
      }
      
      .btn-primary:active {
        transform: translateY(0px) !important;
        box-shadow: 
          0 2px 10px rgba(158, 202, 214, 0.3),
          inset 0 1px 3px rgba(0,0,0,0.2) !important;
      }
      
      /* Luxurious download buttons */
      .btn-default {
        background-color: rgba(42, 31, 53, 0.6) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        color: #9ECAD6 !important;
        border-radius: 6px !important;
        padding: 8px 16px !important;
        transition: all 0.3s ease !important;
        font-family: 'Inter', sans-serif !important;
        font-weight: 300 !important;
      }
      
      .btn-default:hover {
        background-color: rgba(42, 31, 53, 0.9) !important;
        color: #9ECAD6 !important;
        border-color: rgba(158, 202, 214, 0.5) !important;
        box-shadow: 0 4px 12px rgba(158, 202, 214, 0.2) !important;
      }
      
      /* Luxurious tab styling */
      .nav-tabs {
        border-bottom: 1px solid rgba(158, 202, 214, 0.3) !important;
        background-color: transparent !important;
        border-radius: 0 !important;
      }
      
      .nav-tabs > li > a {
        background-color: rgba(42, 31, 53, 0.3) !important;
        color: #9ECAD6 !important;
        border: none !important;
        border-bottom: 2px solid transparent !important;
        border-radius: 0 !important;
        margin-right: 2px !important;
        padding: 12px 20px !important;
        transition: all 0.3s ease !important;
        font-family: 'Inter', sans-serif !important;
        font-weight: 400 !important;
      }
      
      .nav-tabs > li > a:hover {
        background-color: rgba(42, 31, 53, 0.5) !important;
        color: #9ECAD6 !important;
        border-bottom: 2px solid rgba(158, 202, 214, 0.5) !important;
      }
      
      .nav-tabs > li.active > a {
        background-color: rgba(42, 31, 53, 0.7) !important;
        color: #9ECAD6 !important;
        border-bottom: 2px solid rgba(158, 202, 214, 0.8) !important;
        font-weight: 500 !important;
      }
      
      .nav-tabs > li.active > a:hover {
        background-color: rgba(42, 31, 53, 0.8) !important;
        border-bottom: 2px solid rgba(158, 202, 214, 1) !important;
      }
      
      /* Luxurious tab content */
      .tab-content {
        background-color: rgba(26, 21, 32, 0.5) !important;
        border: 1px solid rgba(158, 202, 214, 0.2) !important;
        border-top: none !important;
        border-radius: 0 0 6px 6px !important;
        padding: 25px !important;
      }
      
      /* Luxurious headers */
      h1, h2, h3, h4, h5, h6 {
        font-family: 'Playfair Display', serif !important;
        font-weight: 600 !important;
        letter-spacing: 0.02em !important;
        margin-bottom: 0.75rem !important;
      }
      
      h4, h5 {
        color: #9ECAD6 !important;
        font-weight: 600 !important;
        text-shadow: 0 2px 4px rgba(0, 0, 0, 0.3) !important;
        font-size: 1.3rem !important;
      }
      
      h5 {
        font-size: 1.3rem !important;
        font-weight: 500 !important;
        color: #9ECAD6 !important;
      }
      
      /* Luxurious help text */
      .help-block {
        color: #a89968 !important;
        font-style: italic !important;
        font-size: 15px !important;
        line-height: 1.5 !important;
        font-family: 'Cormorant Garamond', serif !important;
      }
      
      /* Luxurious general text styling */
      p, div, span, label {
        font-family: 'Inter', sans-serif !important;
        font-weight: 300 !important;
        color: #e8e8e8 !important;
      }
      
      /* Luxurious code and monospace elements */
      code, pre, .code {
        font-family: 'Monaco', 'Consolas', monospace !important;
        font-weight: 400 !important;
        font-size: 0.9em !important;
        background: rgba(42, 31, 53, 0.6) !important;
        color: #9ECAD6 !important;
        padding: 4px 8px !important;
        border-radius: 4px !important;
        border: 1px solid rgba(158, 202, 214, 0.2) !important;
        box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.3) !important;
      }
      
      /* Button text */
      .btn {
        font-family: 'Exo 2', sans-serif !important;
        font-weight: 500 !important;
        letter-spacing: 0.01em !important;
        text-transform: none !important;
      }
      
      /* Luxurious input labels */
      .control-label {
        font-family: 'Inter', sans-serif !important;
        font-weight: 400 !important;
        color: #9ECAD6 !important;
        font-size: 16px !important;
        letter-spacing: 0.02em !important;
      }
      
      /* Luxurious horizontal rule */
      hr {
        border-color: rgba(158, 202, 214, 0.2) !important;
        border-width: 1px !important;
        margin: 20px 0 !important;
        box-shadow: 0 1px 0 rgba(158, 202, 214, 0.1) !important;
      }
      
      /* Luxurious data tables */
      .dataTables_wrapper {
        background-color: transparent !important;
        color: #e8e8e8 !important;
      }
      
      .dataTables_wrapper .dataTables_length,
      .dataTables_wrapper .dataTables_filter,
      .dataTables_wrapper .dataTables_info,
      .dataTables_wrapper .dataTables_processing,
      .dataTables_wrapper .dataTables_paginate {
        color: #9ECAD6 !important;
      }
      
      .dataTables_wrapper .dataTables_length select,
      .dataTables_wrapper .dataTables_filter input {
        background-color: rgba(42, 31, 53, 0.6) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        color: #e8e8e8 !important;
        border-radius: 4px !important;
      }
      
      .dataTables_wrapper .dataTables_paginate .paginate_button {
        background-color: rgba(42, 31, 53, 0.5) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        color: #9ECAD6 !important;
        border-radius: 4px !important;
        margin: 0 2px !important;
      }
      
      .dataTables_wrapper .dataTables_paginate .paginate_button:hover {
        background-color: rgba(158, 202, 214, 0.3) !important;
        color: #9ECAD6 !important;
        border-color: rgba(158, 202, 214, 0.6) !important;
      }
      
      .dataTables_wrapper .dataTables_paginate .paginate_button.current {
        background-color: rgba(158, 202, 214, 0.4) !important;
        color: #9ECAD6 !important;
        border-color: rgba(158, 202, 214, 0.7) !important;
      }
      
      /* Luxurious table styling */
      table.dataTable {
        background-color: transparent !important;
        color: #e8e8e8 !important;
      }
      
      table.dataTable thead th {
        background: linear-gradient(135deg, rgba(158, 202, 214, 0.3), rgba(184, 134, 11, 0.3)) !important;
        color: #9ECAD6 !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        font-weight: 500 !important;
        padding: 12px !important;
        font-family: 'Playfair Display', serif !important;
      }
      
      table.dataTable tbody td {
        background-color: rgba(42, 31, 53, 0.3) !important;
        color: #e8e8e8 !important;
        border: 1px solid rgba(158, 202, 214, 0.15) !important;
        padding: 10px !important;
      }
      
      table.dataTable tbody tr:nth-child(even) td {
        background-color: rgba(42, 31, 53, 0.4) !important;
      }
      
      table.dataTable tbody tr:hover td {
        background-color: rgba(158, 202, 214, 0.2) !important;
        color: #ffffff !important;
      }
      
      /* Luxurious regular tables */
      table {
        background-color: transparent !important;
        color: #e8e8e8 !important;
        border-collapse: collapse !important;
        width: 100% !important;
      }
      
      table th {
        background: linear-gradient(135deg, rgba(158, 202, 214, 0.3), rgba(184, 134, 11, 0.3)) !important;
        color: #9ECAD6 !important;
        padding: 12px !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        font-weight: 500 !important;
        font-family: 'Playfair Display', serif !important;
      }
      
      table td {
        background-color: rgba(42, 31, 53, 0.3) !important;
        color: #e8e8e8 !important;
        padding: 10px !important;
        border: 1px solid rgba(158, 202, 214, 0.15) !important;
      }
      
      table tr:nth-child(even) td {
        background-color: rgba(42, 31, 53, 0.4) !important;
      }
      
      /* Luxurious plotly styling */
      .plotly {
        background-color: transparent !important;
      }
      
      /* Luxurious scrollbar styling */
      ::-webkit-scrollbar {
        width: 10px !important;
      }
      
      ::-webkit-scrollbar-track {
        background: rgba(26, 21, 32, 0.5) !important;
      }
      
      ::-webkit-scrollbar-thumb {
        background: rgba(158, 202, 214, 0.5) !important;
        border-radius: 5px !important;
      }
      
      ::-webkit-scrollbar-thumb:hover {
        background: rgba(158, 202, 214, 0.7) !important;
      }
      
      /* Input validation */
      .shiny-input-container {
        color: #e0e0e0 !important;
      }
      
      /* Luxurious select input styling - FIX FOR VISIBILITY */
      select {
        background: rgba(42, 31, 53, 0.9) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        color: #f5f5f5 !important;
        border-radius: 6px !important;
        padding: 10px 14px !important;
        font-family: 'Inter', sans-serif !important;
        font-size: 16px !important;
        font-weight: 300 !important;
        box-shadow: 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 1px 0 rgba(255,255,255,0.03) !important;
        -webkit-appearance: none !important;
        -moz-appearance: none !important;
        appearance: none !important;
      }
      
      select:focus {
        background: rgba(42, 31, 53, 1) !important;
        border-color: rgba(158, 202, 214, 0.6) !important;
        box-shadow: 
          0 0 0 0.2rem rgba(158, 202, 214, 0.1), 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 4px 12px rgba(158, 202, 214, 0.15) !important;
        color: #ffffff !important;
        outline: none !important;
      }
      
      select option {
        background-color: rgba(42, 31, 53, 1) !important;
        color: #f5f5f5 !important;
        padding: 10px 14px !important;
      }
      
      select option:hover {
        background-color: rgba(158, 202, 214, 0.2) !important;
        color: #9ECAD6 !important;
      }
      
      select option:checked {
        background-color: rgba(158, 202, 214, 0.3) !important;
        color: #9ECAD6 !important;
      }
      
       /* Luxurious optimizer dropdown - FIXED VISIBILITY */
         background: rgba(42, 31, 53, 1) !important;
         background-image: none !important;
         border: 1px solid rgba(158, 202, 214, 0.4) !important;
         color: #f5f5f5 !important;
         border-radius: 6px !important;
         padding: 10px 14px !important;
         font-family: 'Inter', sans-serif !important;
         font-size: 16px !important;
         font-weight: 300 !important;
         box-shadow: 
           inset 0 1px 3px rgba(0,0,0,0.3),
           0 1px 0 rgba(255,255,255,0.03) !important;
         -webkit-appearance: menulist !important;
         -moz-appearance: menulist !important;
         appearance: menulist !important;
       }
       
         background: rgba(42, 31, 53, 1) !important;
         background-image: none !important;
         border-color: rgba(158, 202, 214, 0.6) !important;
         box-shadow: 
           0 0 0 0.2rem rgba(158, 202, 214, 0.1), 
           inset 0 1px 3px rgba(0,0,0,0.3),
           0 4px 12px rgba(158, 202, 214, 0.15) !important;
         color: #ffffff !important;
         outline: none !important;
       }
       
         background-color: rgba(42, 31, 53, 1) !important;
         background-image: none !important;
         color: #f5f5f5 !important;
         padding: 10px 14px !important;
         font-family: 'Inter', sans-serif !important;
         border: none !important;
       }
       
         background-color: rgba(158, 202, 214, 0.2) !important;
         background-image: none !important;
         color: #9ECAD6 !important;
       }
       
         background-color: rgba(158, 202, 214, 0.3) !important;
         background-image: none !important;
         color: #9ECAD6 !important;
         font-weight: 400 !important;
       }
       
       /* Luxurious B estimation method dropdown */
        background: rgba(42, 31, 53, 1) !important;
        background-image: none !important;
        border: 1px solid rgba(158, 202, 214, 0.4) !important;
        color: #f5f5f5 !important;
        border-radius: 6px !important;
        padding: 10px 14px !important;
        font-family: 'Inter', sans-serif !important;
        font-size: 16px !important;
        font-weight: 300 !important;
        box-shadow: 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 1px 0 rgba(255,255,255,0.03) !important;
        -webkit-appearance: menulist !important;
        -moz-appearance: menulist !important;
        appearance: menulist !important;
      }
      
        background: rgba(42, 31, 53, 1) !important;
        background-image: none !important;
        border-color: rgba(158, 202, 214, 0.6) !important;
        box-shadow: 
          0 0 0 0.2rem rgba(158, 202, 214, 0.1), 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 4px 12px rgba(158, 202, 214, 0.15) !important;
        color: #ffffff !important;
        outline: none !important;
      }
      
        background-color: rgba(42, 31, 53, 1) !important;
        background-image: none !important;
        color: #f5f5f5 !important;
        padding: 10px 14px !important;
        font-family: 'Inter', sans-serif !important;
        border: none !important;
      }
      
        background-color: rgba(158, 202, 214, 0.2) !important;
        background-image: none !important;
        color: #9ECAD6 !important;
      }
      
        background-color: rgba(158, 202, 214, 0.3) !important;
        background-image: none !important;
        color: #9ECAD6 !important;
        font-weight: 400 !important;
      }
      
      /* Luxurious numeric input styling */
      input[type='number'] {
        background: rgba(42, 31, 53, 0.6) !important;
        border: 1px solid rgba(158, 202, 214, 0.3) !important;
        color: #f5f5f5 !important;
        border-radius: 6px !important;
        padding: 10px 14px !important;
        font-family: 'Inter', sans-serif !important;
        font-size: 16px !important;
        font-weight: 300 !important;
        box-shadow: 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 1px 0 rgba(255,255,255,0.03) !important;
      }
      
      input[type='number']:focus {
        background: rgba(42, 31, 53, 0.8) !important;
        border-color: rgba(158, 202, 214, 0.6) !important;
        box-shadow: 
          0 0 0 0.2rem rgba(158, 202, 214, 0.1), 
          inset 0 1px 3px rgba(0,0,0,0.3),
          0 4px 12px rgba(158, 202, 214, 0.15) !important;
        color: #ffffff !important;
        outline: none !important;
      }
      
       /* Luxurious text input styling */
       input[type='text'] {
         background: rgba(42, 31, 53, 0.6) !important;
         border: 1px solid rgba(158, 202, 214, 0.3) !important;
         color: #f5f5f5 !important;
         border-radius: 6px !important;
         padding: 10px 14px !important;
         font-family: 'Inter', sans-serif !important;
         font-size: 16px !important;
         font-weight: 300 !important;
         box-shadow: 
           inset 0 1px 3px rgba(0,0,0,0.3),
           0 1px 0 rgba(255,255,255,0.03) !important;
       }
       
       input[type='text']:focus {
         background: rgba(42, 31, 53, 0.8) !important;
         border-color: rgba(158, 202, 214, 0.6) !important;
         box-shadow: 
           0 0 0 0.2rem rgba(158, 202, 214, 0.1), 
           inset 0 1px 3px rgba(0,0,0,0.3),
           0 4px 12px rgba(158, 202, 214, 0.15) !important;
         color: #ffffff !important;
         outline: none !important;
       }
       
       /* Luxurious loading spinner styling */
       .shiny-spinner-container {
         background-color: rgba(26, 21, 32, 0.8) !important;
         border: 1px solid rgba(158, 202, 214, 0.3) !important;
         border-radius: 8px !important;
       }
       
       .shiny-spinner-border {
         border-color: rgba(158, 202, 214, 0.3) rgba(158, 202, 214, 0.3) rgba(158, 202, 214, 0.3) rgba(158, 202, 214, 0.8) !important;
         border-width: 3px !important;
       }
       
       .shiny-spinner-text {
         color: #9ECAD6 !important;
         font-family: 'Inter', sans-serif !important;
         font-weight: 300 !important;
         text-shadow: 0 1px 2px rgba(0, 0, 0, 0.5) !important;
       }
       
       /* Custom loading overlay */
       .shiny-spinner-overlay {
         background-color: rgba(10, 10, 10, 0.7) !important;
         backdrop-filter: blur(2px) !important;
       }
       
       /* Luxurious notification styling */
       .shiny-notification {
         background: linear-gradient(135deg, #1a1520, #2a1f35) !important;
         border: 1px solid rgba(158, 202, 214, 0.5) !important;
         border-radius: 8px !important;
         box-shadow: 0 8px 32px rgba(0, 0, 0, 0.6) !important;
         color: #9ECAD6 !important;
         font-family: 'Inter', sans-serif !important;
         font-weight: 400 !important;
       }
       
       .shiny-notification-close {
         color: #9ECAD6 !important;
         opacity: 0.8 !important;
       }
       
       .shiny-notification-close:hover {
         opacity: 1 !important;
         color: #ffffff !important;
       }
       
       /* Success notification */
       .shiny-notification-content .shiny-notification-content {
         background: linear-gradient(135deg, rgba(158, 202, 214, 0.9), rgba(184, 134, 11, 0.9)) !important;
         color: #000000 !important;
       }
     "))
  ),
  div(class = "title-panel",
    div(class = "logo-container",
      div(class = "logo", "⚡"),
      div(class = "logo-text",
        h1("hiddenCOMET", class = "app-title"),
        h2("General CTMC Hidden State Inference", class = "app-subtitle")
      )
    )
  ),
  sidebarLayout(
    sidebarPanel(
      class = "sidebar-panel",
      h4("Process Configuration"),
      textInput("observedStateNames", "Observed state names (auto-detected from B matrix, or specify manually)", value = "state1,state2,state3", placeholder = "e.g., Stage1,Stage2,Stage3"),
      h5("Hidden states: Fixed as E, H, M"),
      hr(),
      h4("Upload inputs"),
      helpText("B: n_observed × 3 CSV emission matrix (any number of observed states × 3 hidden states)."),
      fileInput("fileB", "Emission matrix B (CSV)", accept = c(".csv")),
      helpText("p0: CSV with state names and initial probabilities."),
      fileInput("fileP0", "Initial hidden state distribution p0 (CSV)", accept = c(".csv")),
      helpText("Θ: n_observed × T CSV with timepoints as columns."),
      fileInput("fileTheta", "Observed process Θ (CSV)", accept = c(".csv")),
      hr(),
      numericInput("Npseudo", "Pseudo-count N (for multinomial likelihood)", value = 3000, min = 100, step = 100),
      
      tags$head(
        tags$style(HTML("
      /* Selected value */
      .selectize-control .selectize-input .item { color: #000 !important; }
      .selectize-control .selectize-input input { color: #000 !important; }

      /* Dropdown options */
      .selectize-dropdown .option { color: #000 !important; }
    "))
      ),
      
      numericInput("starts", "Random restarts", value = 30, min = 5, max = 200, step = 5),
      actionButton("runFit", "Submit", class = "btn-primary"),
      width = 4
    ),
    mainPanel(
      class = "main-panel",
      tabsetPanel(
        tabPanel("Inputs",
                 h4("B (emission matrix)"),
                 DTOutput("tblB"),
                 h4("Initial hidden state distribution (p0)"),
                 DTOutput("tblP0"),
                 h4("Observed process Θ (first 10 cols)"),
                 DTOutput("tblThetaHead")
        ),
         tabPanel("Results: Rates & Generator",
                  h4("Estimated generator G"),
                  withSpinner(DTOutput("tblG"), color = "#9ECAD6", type = 4),
                  fluidRow(
                    column(6, h5("Transition rates (per day)"),
                           withSpinner(tableOutput("tblRates"), color = "#9ECAD6", type = 4)),
                    column(6, h5("Derived summaries"),
                           withSpinner(tableOutput("tblSummaries"), color = "#9ECAD6", type = 4))
                  ),
                  downloadButton("downloadG", "Download G (CSV)")
         ),
        tabPanel("Hidden state trajectory p(t)",
                 withSpinner(plotlyOutput("pltEMT", height = "480px"), color = "#9ECAD6", type = 4),
                 downloadButton("downloadPTraj", "Download p(t) (CSV)")
        ),
        tabPanel("Predicted vs Observed Θ",
                 withSpinner(plotlyOutput("pltTheta", height = "520px"), color = "#9ECAD6", type = 4),
                 h5("Prediction residuals (L1 per time)"),
                 withSpinner(tableOutput("tblResiduals"), color = "#9ECAD6", type = 4),
                 downloadButton("downloadThetaPred", "Download Θ_pred (CSV)")
        ),
        tabPanel("Estimate B Matrix",
                 h4("Upload data for B matrix estimation"),
                 fluidRow(
                   column(6,
                   h5("Hidden states trajectory data (P)"),
                   helpText("CSV with 4 columns: timepoints, Epithelial, Hybrid, Mesenchymal"),
                   fileInput("fileEMT", "Hidden states trajectory (CSV)", accept = c(".csv")),
                     numericInput("maxDayEMT", "Maximum day to include", value = 7, min = 1, max = 30)
                   ),
                 column(6,
                   h5("Observed states data (Θ)"),
                   helpText("CSV with timepoints column + any number of observed state columns (e.g., time, state1, state2, state3, state4)"),
                   fileInput("fileCC", "Observed states data (CSV)", accept = c(".csv"))
                 )
                 ),
                 hr(),
                 fluidRow(
                   column(4,
                     selectInput("bMethod", "Estimation method", 
                               choices = c("Least Squares" = "ls", "Maximum Likelihood" = "mle"),
                               selected = "mle")
                   ),
                   column(4,
                     numericInput("bNpseudo", "Pseudo-count N", value = 2000, min = 100, step = 100)
                   ),
                   column(4,
                     numericInput("bStarts", "Random restarts", value = 25, min = 5, max = 100, step = 5)
                   )
                 ),
                  actionButton("runBEstimation", "Estimate B Matrix", class = "btn-primary"),
                  hr(),
                  h4("Estimated B Matrix"),
                  withSpinner(DTOutput("tblBEst"), color = "#9ECAD6", type = 4),
                  br(),
                  withSpinner(plotlyOutput("pltBEst", height = "400px"), color = "#9ECAD6", type = 4),
                  br(),
                  downloadButton("downloadBEst", "Download B Matrix (CSV)")
        ),
        tabPanel("Help / Formats",
                 h4("Expected CSV formats"),
                 tags$ul(
                   tags$li(HTML("<b>B</b>: n_observed × 3 matrix; rownames are observed state names; colnames are the 3 hidden state names. Example for 3×3:
<pre>
,G1,S,G2M
E,0.60,0.30,0.10
H,0.40,0.30,0.30
M,0.20,0.20,0.60
</pre> (If your CSV is transposed, the app detects and corrects.)")),
                   tags$li(HTML("<b>p0</b>: either two columns (State,Value) with hidden state names OR one row with columns matching hidden state names. Values must sum to 1.")),
                   tags$li(HTML("<b>Θ</b>: n_observed × T; rownames are observed state names; columns are numeric timepoints (e.g., 0,0.33,1,3,7). Entries are fractions per column; they will be normalized.")),
      tags$li(HTML("<b>Hidden state trajectory</b>: CSV with 4 columns: timepoints, Epithelial, Hybrid, Mesenchymal")),
      tags$li(HTML("<b>Observed states data</b>: CSV with 4 columns: timepoints, state1, state2, state3 (e.g., time, G1, S, G2M)"))
                 ),
                 h4("Configuration Options"),
                 tags$ul(
                   tags$li("Number of states: Set the number of observed states (2-10). Hidden states are fixed at 3."),
                   tags$li("State names: Customize names for your states (comma-separated)."),
                 ),
                 h4("Notes"),
                 tags$ul(
                   tags$li("The CTMC uses a linear topology where states transition sequentially (E ↔ H ↔ M)."),
                   tags$li("We fit generator by maximizing the multinomial likelihood of Θ given B, p0 via random-restart optimization."),
                   tags$li("Time is taken from Θ column names unless overridden."),
                   tags$li("B matrix estimation uses either Least Squares or Maximum Likelihood methods.")
                 ),
                 br(),
                 p("Please contact annicenajafi27@gmail.com for technical questions.")
      
                 
        )
      )
    )
  )
)

server <- function(input, output, session){
  
  rx_observed_states <- reactive({
    if(!is.null(input$fileB)) {
      tryCatch({
        df <- read.csv(input$fileB$datapath, check.names = FALSE) |> clean_names_safe()
        m <- as.matrix(df)
        if(!is.numeric(m[1,1]) && any(colnames(df) != "")){
          rn <- df[[1]]; df <- df[,-1, drop=FALSE]
          rownames(df) <- rn; m <- as.matrix(df)
        }
        if(!is.null(rownames(m)) && all(rownames(m) != "")) {
          return(rownames(m))
        }
      }, error = function(e) {
      })
    }
    
    states <- trimws(strsplit(input$observedStateNames, ",")[[1]])
    states
  })
  
  rx_hidden_states <- reactive({
    c("E", "H", "M")
  })
  
  iv <- InputValidator$new()
  iv$add_rule("fileB", sv_required(message = "Please upload B (CSV)."))
  iv$add_rule("fileP0", sv_required(message = "Please upload p0 (CSV)."))
  iv$add_rule("fileTheta", sv_required(message = "Please upload Θ (CSV)."))
  iv$enable()
  
  observe({
    cat("Validation status:", iv$is_valid(), "\n")
    if(!iv$is_valid()) {
      cat("Validation errors:\n")
      print(iv$errors)
    }
  })
  
  iv_b <- InputValidator$new()
  iv_b$add_rule("fileEMT", sv_required(message = "Please upload hidden states trajectory data (CSV)."))
  iv_b$add_rule("fileCC", sv_required(message = "Please upload observed states data (CSV)."))
  iv_b$enable()
  
  rx_B <- reactive({
    req(input$fileB)
    df <- read.csv(input$fileB$datapath, check.names = FALSE) |> clean_names_safe()
    m <- as.matrix(df)
    if(!is.numeric(m[1,1]) && any(colnames(df) != "")){
      rn <- df[[1]]; df <- df[,-1, drop=FALSE]
      rownames(df) <- rn; m <- as.matrix(df)
    }
    
     obs_states <- rx_observed_states()
     hidden_states <- rx_hidden_states()
     n_obs <- length(obs_states)
     
     if(nrow(m) != n_obs || ncol(m) != 3){
       stop(paste0("B matrix must be ", n_obs, "×3 (", n_obs, " observed states × 3 hidden states). Found ", nrow(m), "×", ncol(m)))
     }
    
    rownames(m) <- obs_states
    colnames(m) <- hidden_states
    
    for(j in seq_len(ncol(m))){
      colj <- pmax(as.numeric(m[,j]), 0)
      m[,j] <- proj_simplex(colj)
    }
    m
  })
  
  rx_p0 <- reactive({
    req(input$fileP0)
    df <- read.csv(input$fileP0$datapath, check.names = FALSE)
    hidden_states <- rx_hidden_states()
    
    if(all(c("State","Value") %in% names(df))){
      v <- setNames(df$Value, df$State)
    } else {
      if(nrow(df) == 1 && all(hidden_states %in% names(df))){
        v <- setNames(sapply(hidden_states, function(s) df[[s]][1]), hidden_states)
      } else {
        stop(paste0("p0 CSV not recognized. Provide (State,Value) or one row with columns: ", paste(hidden_states, collapse=",")))
      }
    }
    
    missing_states <- setdiff(hidden_states, names(v))
    if(length(missing_states) > 0){
      stop(paste0("Missing states in p0: ", paste(missing_states, collapse=", ")))
    }
    
    v <- pmax(as.numeric(v[hidden_states]), 0)
    v <- proj_simplex(v)
    names(v) <- hidden_states
    v
  })
  
  rx_Theta <- reactive({
    req(input$fileTheta)
    df <- read.csv(input$fileTheta$datapath, check.names = FALSE)
    obs_states <- rx_observed_states()
    n_obs <- length(obs_states)
    
    cat("Theta CSV loaded:", nrow(df), "rows,", ncol(df), "columns\n")
    cat("Observed states:", obs_states, "\n")
    cat("Expected observed states:", n_obs, "\n")
    
    
    time_col <- names(df)[1]
    state_cols <- names(df)[-1]  # Remaining columns are observed states
    
    if(length(state_cols) != n_obs) {
      stop(paste0("CSV should have ", n_obs, " state columns, but found ", length(state_cols), ". Expected states: ", paste(obs_states, collapse = ", ")))
    }
    
    times <- df[[time_col]]
    Theta <- as.matrix(df[, state_cols])
    
    Theta <- t(Theta)
    
    cat("After transposing - Matrix dimensions:", nrow(Theta), "x", ncol(Theta), "\n")
    
    rownames(Theta) <- obs_states
    colnames(Theta) <- as.character(times)
    
    Theta <- apply(Theta, 2, function(x){
      x <- pmax(as.numeric(x), 0); s <- sum(x)
      if(s == 0) rep(1/n_obs, n_obs) else x/s
    })
    rownames(Theta) <- obs_states
    
    Theta
  })
  
  rx_times <- reactive({
    th <- rx_Theta()
    tnames <- colnames(th)
    tvec <- as.numeric(tnames)
    tvec <- sort(unique(tvec))
    validate(need(length(tvec) >= 2, "Need at least two timepoints."))
    tvec
  })
  
  rx_EMT_data <- reactive({
    req(input$fileEMT)
    df <- read.csv(input$fileEMT$datapath, check.names = FALSE)
    
    cat("EMT data loaded:", nrow(df), "rows,", ncol(df), "columns\n")
    cat("EMT column names:", paste(names(df), collapse = ", "), "\n")
    
    if(ncol(df) != 4) {
      stop("EMT trajectory data must have exactly 4 columns: timepoints, Epithelial, Hybrid, Mesenchymal")
    }
    
    original_names <- names(df)
    time_col <- original_names[1]
    emt_cols <- original_names[-1]
    
    expected_cols <- c("Epithelial", "Hybrid", "Mesenchymal")
    names(df)[2:4] <- expected_cols
    emt_cols <- expected_cols
    
    max_day <- input$maxDayEMT
    time_col_name <- names(df)[1]
    if(is.numeric(df[[1]])) {
      df <- df[df[[time_col_name]] <= max_day, ]
    }
    
    times <- df[[1]]
    P <- as.matrix(df[, expected_cols])
    
    cat("After filtering - times length:", length(times), "P dimensions:", dim(P), "\n")
    
    P <- apply(P, 1, function(x) {
      x <- pmax(x, 0)
      if(sum(x) == 0) rep(1/3, 3) else x/sum(x)
    })
    
    cat("After normalization - P dimensions:", dim(P), "times length:", length(times), "\n")
    
    if(ncol(P) == length(times)) {
      colnames(P) <- as.character(times)
      rownames(P) <- expected_cols
    } else {
      stop("Dimension mismatch: P has ", ncol(P), " columns but times has ", length(times), " elements")
    }
    
    P
  })
  
  rx_CC_data <- reactive({
    req(input$fileCC)
    df <- read.csv(input$fileCC$datapath, check.names = FALSE)
    
    cat("CC data loaded:", nrow(df), "rows,", ncol(df), "columns\n")
    cat("CC column names:", paste(names(df), collapse = ", "), "\n")
    
    if(ncol(df) < 3) {
      stop("Observed states data must have at least 3 columns: timepoints + at least 2 observed states")
    }
    
    original_names <- names(df)
    time_col <- original_names[1]  # First column should be timepoints
    state_cols <- original_names[-1]  # Remaining columns are observed states
    
    times <- df[[1]]  # Use first column directly
    Theta <- as.matrix(df[, -1])  # Use all columns except first (timepoints)
    
    cat("CC data - times length:", length(times), "Theta dimensions:", dim(Theta), "\n")
    
    n_states <- ncol(Theta)
    Theta <- apply(Theta, 1, function(x) {
      x <- pmax(x, 0)
      if(sum(x) == 0) rep(1/n_states, n_states) else x/sum(x)
    })
    
    cat("CC data after normalization - Theta dimensions:", dim(Theta), "times length:", length(times), "\n")
    
    if(ncol(Theta) == length(times)) {
      colnames(Theta) <- as.character(times)
      rownames(Theta) <- state_cols
    } else {
      stop("Dimension mismatch: Theta has ", ncol(Theta), " columns but times has ", length(times), " elements")
    }
    
    Theta
  })
  
   b_estimation_state <- eventReactive(input$runBEstimation, {
     isolate({
       req(iv_b$is_valid())
       
      P <- rx_EMT_data()
      Theta <- rx_CC_data()
      
      cat("EMT matrix dimensions:", nrow(P), "x", ncol(P), "\n")
      cat("CC matrix dimensions:", nrow(Theta), "x", ncol(Theta), "\n")
      cat("EMT timepoints:", paste(colnames(P), collapse = ", "), "\n")
      cat("CC timepoints:", paste(colnames(Theta), collapse = ", "), "\n")
      
      p_times <- colnames(P)
      theta_times <- colnames(Theta)
      
      if(!setequal(p_times, theta_times)) {
        stop("Timepoints in EMT trajectory and observed states data must match exactly. 
             EMT times: ", paste(p_times, collapse = ", "), 
             " | Observed times: ", paste(theta_times, collapse = ", "))
      }
      
      common_times <- as.character(sort(as.numeric(p_times)))
      validate(need(length(common_times) >= 2, "Need at least two timepoints."))
       
       P <- P[, common_times, drop = FALSE]
       Theta <- Theta[, common_times, drop = FALSE]
       
       showNotification("Estimating B matrix...", duration = NULL, id = "b_estimation")
       
      if(input$bMethod == "ls") {
        cat("Using Least Squares method\n")
        B_hat <- estimate_B_ls(Theta, P)
      } else {
        cat("Using Maximum Likelihood method\n")
        B_hat <- estimate_B_mle(Theta, P, observed_states = rownames(Theta), hidden_states = c("E", "H", "M"), N = input$bNpseudo, starts = input$bStarts)
      }
      
      cat("Final B matrix dimensions:", dim(B_hat), "\n")
       
       removeNotification("b_estimation")
       showNotification("B matrix estimation completed!", type = "message", duration = 3)
       
       list(B = B_hat, method = input$bMethod, times = common_times)
     })
   }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$tblB <- renderDT({
    req(rx_B())
    datatable(as.data.frame(rx_B()) |> rownames_to_column("ObservedState"),
              options = list(
                dom = 't',
                pageLength = 10,
                scrollX = TRUE,
                initComplete = JS("function(settings, json) { $(this.api().table().container()).css({'background-color': '#2d2d2d', 'color': '#e0e0e0'}); }")
              ), 
              rownames = FALSE,
              class = "display compact hover")
  })
  
  output$tblP0 <- renderDT({
    req(rx_p0())
    p0 <- rx_p0()
    hidden_states <- names(p0)
    datatable(data.frame(State = hidden_states, Value = rnd(p0,3)),
              options = list(
                dom = 't',
                initComplete = JS("function(settings, json) { $(this.api().table().container()).css({'background-color': '#2d2d2d', 'color': '#e0e0e0'}); }")
              ), 
              rownames = FALSE,
              class = "display compact hover")
  })
  
  output$tblThetaHead <- renderDT({
    req(rx_Theta())
    th <- rx_Theta()
    if(ncol(th) > 10) th <- th[, 1:10, drop=FALSE]
    datatable(as.data.frame(th) |> rownames_to_column("ObservedState"),
              options = list(
                scrollX = TRUE, 
                dom = 't',
                initComplete = JS("function(settings, json) { $(this.api().table().container()).css({'background-color': '#2d2d2d', 'color': '#e0e0e0'}); }")
              ), 
              rownames = FALSE,
              class = "display compact hover")
  })
  
  fitting_results <- reactiveValues(data = NULL)
  
  observeEvent(input$runFit, {
    cat("fit_state eventReactive triggered!\n")
    isolate({
      tryCatch({
        cat("Starting fitting process...\n")
        if(!iv$is_valid()) {
          cat("Validation failed, stopping\n")
          return(NULL)
        }
        cat("Validation passed, reading data...\n")
       B <- rx_B()
       p0 <- rx_p0()
       Theta <- rx_Theta()
       times <- rx_times()
       hidden_states <- rx_hidden_states()
       
       cat("Data loaded successfully\n")
       cat("B dimensions:", dim(B), "\n")
       cat("Theta dimensions:", dim(Theta), "\n")
       cat("Times:", times, "\n")
       
       Theta <- Theta[, match(times, as.numeric(colnames(Theta))), drop=FALSE]
       B <- proj_columns_simplex(B)
       
       showNotification("Fitting generator matrix...", duration = NULL, id = "fitting")
       
       cat("Starting optimization...\n")
       res <- fit_G_given_p0(
         Theta = Theta,
         B = B,
         times = times,
         p0 = p0,
         states = hidden_states,
         starts = input$starts,
         N_pseudo = input$Npseudo,
         method = "Nelder-Mead"
       )
       cat("Optimization completed successfully\n")
       
      removeNotification("fitting")
      showNotification("Fitting completed!", type = "message", duration = 3)

      fitting_results$data <- list(B = B, p0 = res$p0, G = res$G, P = res$P, Theta = Theta, times = times)

      updateTabsetPanel(session, "tabs", selected = "Results: Rates & Generator")

      }, error = function(e) {
        cat("Error in fitting:", e$message, "\n")
        removeNotification("fitting")
        showNotification(paste("Error:", e$message), type = "error", duration = 5)
        fitting_results$data <- NULL
      })
    })
  })
  
  fit_state <- reactive({
    fitting_results$data
  })
  
  output$tblG <- renderDT({
    fs <- fit_state(); req(fs)
    datatable(as.data.frame(fs$G) |> rownames_to_column("from"),
              options = list(
                dom = 't',
                initComplete = JS("function(settings, json) { $(this.api().table().container()).css({'background-color': '#2d2d2d', 'color': '#e0e0e0'}); }")
              ), 
              rownames = FALSE,
              class = "display compact hover")
  })
  
  output$tblRates <- renderTable({
    fs <- fit_state(); req(fs)
    G <- fs$G
    hidden_states <- rx_hidden_states()
    
    rates <- data.frame()
    for(i in seq_len(nrow(G))){
      for(j in seq_len(ncol(G))){
        if(i != j && G[i,j] > 0){
          from_state <- rownames(G)[i]
          to_state <- colnames(G)[j]
          rate <- G[i,j]
          rates <- rbind(rates, data.frame(
            Parameter = paste0(from_state, "→", to_state),
            Rate_per_day = round(rate, 4)
          ))
        }
      }
    }
    rates
  })
  
  output$tblSummaries <- renderTable({
    fs <- fit_state(); req(fs)
    G <- fs$G
    hidden_states <- rx_hidden_states()
    
    leave <- -diag(G)
    halft <- log(2) / leave
    names(leave) <- hidden_states
    names(halft) <- hidden_states
    
    A <- t(G); A[ , ncol(A)] <- 1
    b <- c(rep(0, nrow(G)-1), 1)
    pi <- MASS::ginv(A) %*% b
    pi <- as.numeric(pi); pi <- proj_simplex(pi)
    names(pi) <- hidden_states
    
    data.frame(
      Quantity = c(paste("Mean leaving rate", hidden_states),
                   paste("Half-time", hidden_states, "(days)"),
                   paste("Stationary", hidden_states)),
      Value = c(leave, halft, pi) |> round(4)
    )
  })
  
  output$pltEMT <- renderPlotly({
    fs <- fit_state(); req(fs)
    P <- fs$P; times <- fs$times
    hidden_states <- rx_hidden_states()
    
    df <- as_tibble(P, rownames = "State") |>
      pivot_longer(-State, names_to = "time", values_to = "frac") |>
      mutate(time = as.numeric(time),
             State = factor(State, levels = hidden_states))
    
    n_states <- length(hidden_states)
    colors <- c('#C5B0CD','#FFB8E0','#2F5755')
    names(colors) <- hidden_states
    
    p <- ggplot(df, aes(time, frac, color = State)) +
      geom_line(size = 1.2) + geom_point(size = 2, stroke=1.5, shape=8) +
      scale_color_manual(values = colors) +
      labs(title = "Inferred hidden state trajectory p(t)",
           x = "Time (days)", y = "Fraction") +
          theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, color = "#B91C5C"),
        plot.background = element_rect(fill = "#2d2d2d", color = NA),
        panel.background = element_rect(fill = "#2d2d2d", color = NA),
        panel.grid.major = element_line(color = "#4d4d4d", size = 0.5),
        panel.grid.minor = element_line(color = "#3d3d3d", size = 0.3),
        axis.text = element_text(color = "#e0e0e0"),
        axis.title = element_text(color = "#e0e0e0"),
        legend.background = element_rect(fill = "#2d2d2d", color = "#8B1538"),
        legend.text = element_text(color = "#e0e0e0"),
        legend.title = element_text(color = "#B91C5C")
      )
    ggplotly(p)
  })
  
  output$pltTheta <- renderPlotly({
    fs <- fit_state(); req(fs)
    B <- fs$B; P <- fs$P; Theta <- fs$Theta
    observed_states <- rx_observed_states()
    Th_pred <- B %*% P
    
    df_obs <- as_tibble(Theta, rownames = "Phase") |>
      pivot_longer(-Phase, names_to = "time", values_to = "value") |>
      mutate(time = as.numeric(time), Kind = "Observed")
    df_pred <- as_tibble(Th_pred, rownames = "Phase") |>
      pivot_longer(-Phase, names_to = "time", values_to = "value") |>
      mutate(time = as.numeric(time), Kind = "Predicted")
    df <- bind_rows(df_obs, df_pred)
    
    p <- ggplot(df, aes(time, value, color = Kind)) +
      geom_line(size = 1.1) + geom_point(size = 2) +
      facet_wrap(~Phase, nrow = 1) +
      scale_color_manual(values = c(Observed="#B91C5C", Predicted="#8B1538")) +
      labs(title = "Observed process Θ: observed vs. model prediction",
           x = "Time (days)", y = "Fraction") +
          theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, color = "#B91C5C"),
        plot.background = element_rect(fill = "#2d2d2d", color = NA),
        panel.background = element_rect(fill = "#2d2d2d", color = NA),
        panel.grid.major = element_line(color = "#4d4d4d", size = 0.5),
        panel.grid.minor = element_line(color = "#3d3d3d", size = 0.3),
        axis.text = element_text(color = "#e0e0e0"),
        axis.title = element_text(color = "#e0e0e0"),
        legend.background = element_rect(fill = "#2d2d2d", color = "#8B1538"),
        legend.text = element_text(color = "#e0e0e0"),
        legend.title = element_text(color = "#B91C5C"),
        strip.background = element_rect(fill = "#8B1538", color = "#B91C5C"),
        strip.text = element_text(color = "#ffffff", face = "bold")
      )
    ggplotly(p)
  })
  
  output$tblResiduals <- renderTable({
    fs <- fit_state(); req(fs)
    B <- fs$B; P <- fs$P; Theta <- fs$Theta
    Th_pred <- B %*% P
    times <- as.numeric(colnames(Theta))
    L1 <- colSums(abs(Theta - Th_pred))
    data.frame(Time = times, L1_residual = round(L1, 4))
  })
  
  output$downloadG <- downloadHandler(
    filename = function(){ "G_estimated.csv" },
    content = function(file){
      fs <- fit_state(); req(fs)
      write.csv(fs$G, file)
    }
  )
  
  output$downloadPTraj <- downloadHandler(
    filename = function(){ "EMT_trajectory_p_t.csv" },
    content = function(file){
      fs <- fit_state(); req(fs)
      write.csv(fs$P, file)
    }
  )
  
  output$downloadThetaPred <- downloadHandler(
    filename = function(){ "Theta_predicted.csv" },
    content = function(file){
      fs <- fit_state(); req(fs)
      Th_pred <- fs$B %*% fs$P
      write.csv(Th_pred, file)
    }
  )
  
  output$tblBEst <- renderDT({
    bes <- b_estimation_state(); req(bes)
    datatable(as.data.frame(bes$B) |> rownames_to_column("CyclePhase"),
              options = list(
                dom = 't',
                initComplete = JS("function(settings, json) { $(this.api().table().container()).css({'background-color': '#2d2d2d', 'color': '#e0e0e0'}); }")
              ), 
              rownames = FALSE,
              class = "display compact hover",
              caption = paste("Estimated B matrix using", ifelse(bes$method == "ls", "Least Squares", "Maximum Likelihood"), "method"))
  })
  
  output$pltBEst <- renderPlotly({
    bes <- b_estimation_state(); req(bes)
    B <- bes$B
    
    cat("B matrix in plot:", dim(B), "\n")
    cat("B matrix rownames:", rownames(B), "\n")
    cat("B matrix colnames:", colnames(B), "\n")
    
    df <- as_tibble(B, rownames = "Observed") |>
      pivot_longer(-Observed, names_to = "Hidden", values_to = "prob") |>
      mutate(
        Observed = factor(Observed, levels = rownames(B)),
        Hidden = factor(Hidden, levels = colnames(B))
      )
    
    cat("Plot data dimensions:", nrow(df), "rows\n")
    cat("Unique observed states:", unique(df$Observed), "\n")
    cat("Unique hidden states:", unique(df$Hidden), "\n")
    
    pink_purple_strong <- c("#FCE4EC", "#F48FB1", "#CE5ABD", "#8E24AA", "#4A148C", "#1A0033")
    
    p <- ggplot(df, aes(x = Hidden, y = Observed, fill = prob, text = paste0("Observed: ", Observed, "<br>Hidden: ", Hidden, "<br>Probability: ", round(prob, 3)))) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", prob)), color = "white", size = 4) +
      scale_fill_gradientn(
        colors = pink_purple_strong,
        limits = c(0, 1),
        values = scales::rescale(c(0, 0.15, 0.35, 0.6, 0.85, 1)),
        name = "Probability"
      ) +
      labs(
        title = paste("Estimated B Matrix (", ifelse(bes$method == "ls", "Least Squares", "Maximum Likelihood"), ")"),
        x = "Hidden State", 
        y = "Observed State"
      ) +
          theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, color = "#B91C5C"),
        plot.background = element_rect(fill = "#2d2d2d", color = NA),
        panel.background = element_rect(fill = "#2d2d2d", color = NA),
        axis.text = element_text(color = "#e0e0e0"),
        axis.title = element_text(color = "#e0e0e0"),
        legend.background = element_rect(fill = "#2d2d2d", color = "#8B1538"),
        legend.text = element_text(color = "#e0e0e0"),
        legend.title = element_text(color = "#B91C5C"),
        panel.grid = element_blank()
      )
    
    ggplotly(p, tooltip = "text") |>
      layout(
        plot_bgcolor = "#2d2d2d",
        paper_bgcolor = "#2d2d2d",
        font = list(color = "#e0e0e0")
      )
  })
  
  output$downloadBEst <- downloadHandler(
    filename = function(){ 
      paste0("B_matrix_estimated_", input$bMethod, "_", input$cellLineName, ".csv") 
    },
    content = function(file){
      bes <- b_estimation_state(); req(bes)
      write.csv(bes$B, file)
    }
  )
}

shinyApp(ui, server)
