score_statistic <- function(Un,
                            Kn,
                            indices_H0) {
  # Subvector of the score function
  U_psi <- Un[indices_H0]

  # [K^-1]_{psi,psi}: invert the FULL Fisher information first, then take the
  # submatrix at indices_H0 (the Schur complement of the nuisance-parameter
  # block) - NOT the inverse of the submatrix, which coincides with this only
  # when the tested parameters are Fisher-orthogonal to the nuisance ones.
  # Unlike the Wald statistic, LM = U_psi' [K^-1]_{psi,psi} U_psi uses this
  # block directly - no further inversion.
  R <- diag(length(Un))[indices_H0, , drop = FALSE]
  K_inv_psi_psi <- R %*% tryCatch(
    solve(Kn, t(R)),
    error = function(e) MASS::ginv(Kn) %*% t(R)
  )

  # Score statistic and p-value
  S <- as.numeric(t(U_psi) %*% K_inv_psi_psi %*% U_psi)
  list(statistic = S, p_value = 1 - pchisq(S, df = length(indices_H0)))
}