wald_statistic <- function(theta_hat,  # Estimated parameters
                           Kn,         # Fisher information matrix
                           theta0,     # Hypothesized values under H0
                           indices_H0  # Indices of parameters under H0
                           ) {
  # Matrix for the hypothesis test
  R <- diag(length(theta_hat))[indices_H0, , drop = FALSE]

  # Difference between estimated and hypothesized parameters
  delta <- R %*% theta_hat - theta0

  # Variance of the linear combination R*theta
  K_psi_psi <- R %*% tryCatch(
    solve(Kn, t(R)),
    error = function(e) MASS::ginv(Kn) %*% t(R)
  )

  K_psi_inv <- tryCatch(
    solve(K_psi_psi),
    error = function(e) MASS::ginv(K_psi_psi)
  )

  # Wald statistic and p-value
  W <- as.numeric(t(delta) %*% K_psi_inv %*% delta)
  list(statistic = W, p_value = 1 - pchisq(W, df = length(indices_H0)))
}