score_statistic <- function(Un,
                            Kn,
                            indices_H0) {
  # Subvector of the score function
  U_psi <- Un[indices_H0]

  # Inverse of the submatrix of Fisher information
  K_psi <- Kn[indices_H0, indices_H0, drop = FALSE]
  K_psi_inv <- tryCatch(solve(K_psi), error = function(e) MASS::ginv(K_psi))

  # Score statistic and p-value
  S <- as.numeric(t(U_psi) %*% K_psi_inv %*% U_psi)
  list(statistic = S, p_value = 1 - pchisq(S, df = length(indices_H0)))
}