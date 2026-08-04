robust_wald_statistic <- function(theta_hat,
                                  Hn,
                                  Jn,
                                  theta0,
                                  indices_H0) {
  # Matrix for the hypothesis test
  R <- diag(length(theta_hat))[indices_H0, , drop = FALSE]

  # Difference between estimated and hypothesized parameters
  delta <- R %*% theta_hat - theta0

  # Sandwich covariance matrix restricted to the subvector
  Hn_inv <- tryCatch(solve(Hn), error = function(e) MASS::ginv(Hn))
  Sigmaq <- R %*% Hn_inv %*% Jn %*% Hn_inv %*% t(R)
  Sigmaq_inv <- tryCatch(solve(Sigmaq), error = function(e) MASS::ginv(Sigmaq))

  # Robust Wald statistic and p-value
  W <- as.numeric(t(delta) %*% Sigmaq_inv %*% delta)
  list(statistic = W, p_value = 1 - pchisq(W, df = length(indices_H0)))
}