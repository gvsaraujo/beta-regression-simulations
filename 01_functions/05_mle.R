#**** Estimating the parameters via MLE****
MLE_BETA <- function(y, X, Z, linkmu="logit", linkphi="log", method = "BFGS",
                     beta_fix = NULL, idx_beta = NULL,  
                     gamma_fix = NULL, idx_gamma = NULL) {
  
  
  # y: response variable (numeric vector)
  # X: design matrix for beta (numeric matrix)
  # Z: design matrix for gamma (numeric matrix)
  # linkmu: link function for mu (character, default "logit")
  # linkphi: link function for phi (character, default "log")
  # method: Method for the optimization
  # beta_fix: fixed values for beta parameters (numeric vector)
  # idx_beta: indices of fixed beta parameters (numeric vector)
  # gamma_fix: fixed values for gamma parameters (numeric vector)
  # idx_gamma: indices of fixed gamma parameters (numeric vector)
  
  #-------->Define some quantities<------------#
  n <- length(y) # Number of observations
  p_total <- ncol(X) # Total number of beta parameters
  q_total <- ncol(Z) # Total number of gamma parameters
  
  # Free values
  X_free <- if (is.null(idx_beta)) X else X[ , -idx_beta, drop = FALSE] # Free design matrix for beta
  Z_free <- if (is.null(idx_gamma)) Z else Z[ , -idx_gamma, drop = FALSE] # Free design matrix for gamma
  
  #-------->Define free parameter indices<------------#
  if (!is.null(idx_beta)) {
    idx_beta_free <- setdiff(1:p_total, idx_beta)
  } else {
    idx_beta_free <- 1:p_total
  }
  if (!is.null(idx_gamma)) {
    idx_gamma_free <- setdiff(1:q_total, idx_gamma)
  } else {
    idx_gamma_free <- 1:q_total
  }
  
  #-------->Free parameter counts<------------#
  p_free <- length(idx_beta_free) # Number of free beta parameters
  q_free <- length(idx_gamma_free) # Number of free gamma parameters
  total_par <- p_free + q_free # Total number of parameters to be estimated
  
  #-------->Start estimation process<------------#
  
  if (total_par == 0) { # If no parameters are free, we can directly compute the results.
    
    theta_hat_full <- c(beta_fix, gamma_fix) # Combine fixed parameters
    
    results <- list(theta_hat_full = theta_hat_full, # Full parameter estimates
                    score = Score_LSMLE(theta_hat_full, y=y, X=X, Z=Z, alpha_const = 0, linkmu=linkmu, linkphi=linkphi), # Score vector
                    fisher_info = fisher_info_beta_regression(y, X, Z, theta_hat_full) # Fisher information matrix
    )
    return(results) # Return results
  }
  
  
  #*** Log-likelihood function for optimization ***
  loglik <- function(theta_free) {
    
    beta_full <- numeric(p_total)
    gamma_full <- numeric(q_total)
    
    if (!is.null(beta_fix)) beta_full[idx_beta] <- beta_fix
    if (p_free > 0) beta_full[idx_beta_free] <- theta_free[1:p_free]
    
    if (!is.null(gamma_fix)) gamma_full[idx_gamma] <- gamma_fix
    if (q_free > 0) gamma_full[idx_gamma_free] <- theta_free[(p_free + 1):(p_free + q_free)]
    
    eta_mu <- as.vector(X %*% beta_full)
    eta_phi <- as.vector(Z %*% gamma_full)
    mu <- mean_link(beta = beta_full, y = y, X = X, linkmu = linkmu)$mu
    phi <- precision_link(gama = gamma_full, Z = Z, linkphi = linkphi)$phi
    
    a <- mu * phi
    b <- (1 - mu) * phi
    
    tryCatch({
      ll_components <- lgamma(phi) - lgamma(a) - lgamma(b) +
        (a - 1) * log(y) + (b - 1) * log(1 - y)
      
      ll <- sum(ll_components)
      
      return(ll)
    }, error = function(e) {
      stop("Error in log-likelihood calculation: ", e$message)
    })
  }
  
  
  init_values <- initial_values_beta_model(y, X, Z, linkmu = linkmu, 
                                           idx_beta = idx_beta, idx_gamma = idx_gamma)
  # Optimization
  fit <- optim(par = init_values, fn = loglik,
               method = method, control = list(fnscale = -1,maxit = 5000, reltol = 1e-8))
  
  # Estimated parameters
  theta_hat <- fit$par
  theta_hat_full <- numeric(p_total + q_total)
  
  if (!is.null(beta_fix)) theta_hat_full[idx_beta] <- beta_fix
  if (p_free > 0) theta_hat_full[idx_beta_free] <- theta_hat[1:p_free]
  
  if (!is.null(gamma_fix)) theta_hat_full[p_total + idx_gamma] <- gamma_fix
  if (q_free > 0) theta_hat_full[p_total + idx_gamma_free] <- theta_hat[(p_free + 1):(p_free + q_free)]

  U_hat_full <- Score_LSMLE(theta_hat_full, y=y, X=X, Z=Z, alpha_const = 0, linkmu=linkmu, linkphi=linkphi) # Score vector
  I_hat_full <- fisher_info_beta_regression(y, X, Z, theta_hat_full)

  list(
    theta_hat = theta_hat,
    #theta_hat_full = theta_hat_full,
    theta_hat_full = theta_hat_full,
    score = U_hat_full,
    fisher_info = I_hat_full)
}