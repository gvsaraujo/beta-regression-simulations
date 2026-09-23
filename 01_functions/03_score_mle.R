#***Score vector for MLE****
score_beta_regression <- function(y, X, Z, theta_hat, linkmu = "logit", linkphi = "log") {
  p <- ncol(X)
  q <- ncol(Z)
  beta <- theta_hat[1:p]
  gamma <- theta_hat[(p + 1):(p + q)]
  mean_info <- mean_link(beta = beta, y = y, X = X, linkmu = linkmu)
  mu <- mean_info$mu
  T_1 <- mean_info$T_1
  precision_info <- precision_link(gama = gamma, Z = Z, linkphi = linkphi)
  phi <- precision_info$phi
  T_2 <- precision_info$T_2
  a <- mu * phi
  b <- (1 - mu) * phi
  log_y <- log(y)
  log_1_y <- log(1 - y)
  dig_mu_phi <- digamma(a)
  dig_1mu_phi <- digamma(b)
  dig_phi <- digamma(phi)

  score_mu <- phi * (log_y - log_1_y - dig_mu_phi + dig_1mu_phi)
  grad_beta <- as.vector(crossprod(X, score_mu * T_1))

  score_phi <- mu * (log_y - dig_mu_phi) + (1 - mu) * (log_1_y - dig_1mu_phi) + dig_phi
  grad_gamma <- as.vector(crossprod(Z, score_phi * T_2))

  return(c(grad_beta, grad_gamma))
}