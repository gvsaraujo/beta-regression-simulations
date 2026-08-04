robust_score_statistic <- function(Un,
                                   Hn,
                                   Jn,
                                   indices_H0) {
  # Subvector of the estimating function
  U_psi <- Un[indices_H0]

  # Block of H_inv corresponding to psi
  Hn_inv    <- tryCatch(solve(Hn), error = function(e) MASS::ginv(Hn))
  H_psi_psi <- Hn_inv[indices_H0, indices_H0, drop = FALSE]

  # Block of (Hn %*% Jn_inv %*% Hn) corresponding to psi
  Jn_inv <- tryCatch(solve(Jn),  error = function(e) MASS::ginv(Jn))
  M <- tryCatch(solve(Hn %*% Jn_inv %*% Hn), error = function(e) MASS::ginv(Hn %*% Jn_inv %*% Hn))
  M_psi_inv <- tryCatch(solve(M[indices_H0, indices_H0, drop = FALSE]),
                        error = function(e) MASS::ginv(M[indices_H0, indices_H0, drop = FALSE]))

  # Robust score statistic and p-value
  S <- as.numeric(t(U_psi) %*% H_psi_psi %*% M_psi_inv %*% H_psi_psi %*% U_psi)
  list(statistic = S, p_value = 1 - pchisq(S, df = length(indices_H0)))
}