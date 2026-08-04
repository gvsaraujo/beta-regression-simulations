gradient_statistic <- function(Un,
                               theta_hat,
                               theta0,
                               indices_H0) {
  # Difference between estimated and hypothesized parameters
  R     <- diag(length(Un))[indices_H0, , drop = FALSE]
  delta <- R %*% theta_hat - theta0

  # Gradient statistic and p-value
  G <- as.numeric(Un[indices_H0] %*% delta)
  list(statistic = G, p_value = 1 - pchisq(G, df = length(indices_H0)))
}