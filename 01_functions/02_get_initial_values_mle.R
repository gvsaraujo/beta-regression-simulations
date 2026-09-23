#****Initial values for Beta Regression (MLE) - Ferrari-Cribari Neto***
initial_values_beta_model <- function(y, X, Z, linkmu = "logit",
                                      idx_beta = NULL, idx_gamma = NULL) {
  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  
  # Apply the mean_link
  link_info <- mean_link(beta = rep(0, p), y = y, X = X, linkmu = linkmu)
  yg1 <- link_info$yg1
  
  # Initial values using linear regression
  beta_init <- solve(t(X) %*% X) %*% t(X) %*% yg1
  
  # Calculating the values of mu_hat and eta_hat
  eta_hat <- as.vector(X %*% beta_init)
  mu_hat <- switch(linkmu,
                   "logit"   = exp(eta_hat) / (1 + exp(eta_hat)),
                   "probit"  = pnorm(eta_hat),
                   "cloglog" = 1 - exp(-exp(eta_hat)),
                   "log"     = exp(eta_hat),
                   "loglog"  = exp(-exp(-eta_hat)),
                   "cauchit" = atan(eta_hat) / pi + 0.5)
  
  # Residuals
  residuals <- yg1 - eta_hat
  
  # Derivative of the inverse link evaluated at mu_hat
  g_prime <- link_info$T_1
  
  # Variance estimator
  sigma2 <- (residuals^2) / ((n - p) * (g_prime^2))
  
  # Calculating phi
  phi_vec <- mu_hat * (1 - mu_hat) / sigma2
  phi_scalar <- max(mean(phi_vec) - 1, 1)  # ensure phi > 0
  
  gamma_init <- c(log(phi_scalar), rep(0, q - 1))
  
  # beta free index
  beta_free_idx <- setdiff(seq_along(beta_init), idx_beta)
  initi_beta_free <- beta_init[beta_free_idx]
  
  # gamma free index
  gamma_free_idx <- setdiff(seq_along(gamma_init), idx_gamma)
  initi_gamma_free <- gamma_init[gamma_free_idx]
  
  # Getting initial values for optimization
  initial_vals <- c(initi_beta_free, initi_gamma_free)
  
  return(initial_vals)
}