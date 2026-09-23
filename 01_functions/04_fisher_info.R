#****Fisher information matrix for MLE****
fisher_info_beta_regression <- function(y, X, Z, theta_hat, linkmu = "logit", linkphi = "log") {
  p <- ncol(X)
  q <- ncol(Z)
  beta <- theta_hat[1:p]
  gamma <- theta_hat[(p+1):(p+q)]
  mean_info <- mean_link(beta, y, X, linkmu)
  precision_info <- precision_link(gamma, Z, linkphi)
  mu <- mean_info$mu
  phi <- precision_info$phi
  dmu_deta <- mean_info$T_1
  dphi_deta <- precision_info$T_2
  psi1_mu_phi <- trigamma(mu * phi)
  psi1_1mu_phi <- trigamma((1 - mu) * phi)
  psi1_phi <- trigamma(phi)

  wi_bb <- phi^2 * (psi1_mu_phi + psi1_1mu_phi) * dmu_deta^2
  wi_bg <- phi * (mu * psi1_mu_phi - (1 - mu) * psi1_1mu_phi) * dphi_deta * dmu_deta
  wi_gg <- (mu^2 * psi1_mu_phi + (1 - mu)^2 * psi1_1mu_phi - psi1_phi) * dphi_deta^2

  I_bb <- crossprod(X, X * wi_bb)
  I_bg <- crossprod(X, Z * wi_bg)
  I_gg <- crossprod(Z, Z * wi_gg)

  I <- rbind(cbind(I_bb, I_bg), cbind(t(I_bg), I_gg))
  return(I)
}