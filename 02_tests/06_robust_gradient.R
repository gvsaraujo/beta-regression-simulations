robust_gradient_statistic <- function(Un,
                                      theta_hat,
                                      theta0,
                                      Hn,
                                      Jn,
                                      indices_H0) {
  # Subvector of the estimating function and parameter difference
  U_psi <- Un[indices_H0]
  delta <- theta_hat[indices_H0] - theta0

  # Block of Jn_inv %*% Hn corresponding to psi
  Jn_inv <- tryCatch(solve(Jn), error = function(e) MASS::ginv(Jn))
  M_psi_psi <- (Jn_inv %*% Hn)[indices_H0, indices_H0, drop = FALSE]

  # Robust gradient statistic and p-value
  T_stat <- as.numeric(-t(U_psi) %*% M_psi_psi %*% delta)
  list(statistic = T_stat, p_value = 1 - pchisq(T_stat, df = length(indices_H0)))
}