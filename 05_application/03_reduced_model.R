# Reduced beta regression model for the HIC application (Urb dropped from
# the precision submodel, since Table 4.2 shows it is not significant for
# any estimator). Reproduces Table 4.3 of the dissertation.
#
# Model (all 80 observations, including the influential Observation #1):
#   logit(mu_i)  = beta1 + beta2*Urb_i + beta3*GDP_i
#   log(phi_i)   = gamma1 + gamma2*GDP_i

if (!dir.exists("../01_functions") || !dir.exists("../02_tests")) {
  stop("Run with the working directory set to 05_application/")
}
invisible(lapply(list.files("../01_functions", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("../02_tests", pattern = "\\.R$", full.names = TRUE), source))

d <- read.table("HICdata.txt", header = TRUE)
y <- d$HIC
X <- cbind(1, d$Urb, d$GDP)
Z <- cbind(1, d$GDP)
param_names <- c("beta1 (Int.)", "beta2 (Urb)", "beta3 (GDP)",
                 "gamma1 (Int.)", "gamma2 (GDP)")

fit_mle    <- MLE_BETA(y, X, Z)
fit_lsmle  <- LSMLE_BETA(y, X, Z)
fit_lmdpde <- LMDPDE_BETA(y, X, Z)

se_mle <- sqrt(diag(solve(fit_mle$fisher_info)))
cov_lsmle  <- solve(fit_lsmle$Vq$H)  %*% fit_lsmle$Vq$J  %*% solve(fit_lsmle$Vq$H)
cov_lmdpde <- solve(fit_lmdpde$Vq$H) %*% fit_lmdpde$Vq$J %*% solve(fit_lmdpde$Vq$H)
se_lsmle  <- sqrt(diag(cov_lsmle))
se_lmdpde <- sqrt(diag(cov_lmdpde))

cat("=== Table 4.3 (reduced model): coefficient estimates (SE) ===\n")
print(data.frame(Parameter = param_names,
                 MLE = round(fit_mle$theta_hat_full, 4), MLE_SE = round(se_mle, 4),
                 LSMLE = round(fit_lsmle$theta_hat, 4), LSMLE_SE = round(se_lsmle, 4),
                 LMDPDE = round(fit_lmdpde$theta_hat, 4), LMDPDE_SE = round(se_lmdpde, 4)),
     row.names = FALSE)
cat("alpha_const: LSMLE =", fit_lsmle$alpha_const, " LMDPDE =", fit_lmdpde$alpha_const, "\n\n")

test_rows <- list()
for (i in seq_along(param_names)) {
  is_beta <- i <= 3
  idx <- if (is_beta) i else i - 3

  if (is_beta) {
    fit_mle_r    <- MLE_BETA(y, X, Z, beta_fix = 0, idx_beta = idx)
    fit_lsmle_r  <- LSMLE_BETA(y, X, Z, beta_fix = 0, idx_beta = idx,
                               qoptimal = FALSE, q0 = 1 - fit_lsmle$alpha_const)
    fit_lmdpde_r <- LMDPDE_BETA(y, X, Z, beta_fix = 0, idx_beta = idx,
                                qoptimal = FALSE, q0 = 1 - fit_lmdpde$alpha_const)
  } else {
    fit_mle_r    <- MLE_BETA(y, X, Z, gamma_fix = 0, idx_gamma = idx)
    fit_lsmle_r  <- LSMLE_BETA(y, X, Z, gamma_fix = 0, idx_gamma = idx,
                               qoptimal = FALSE, q0 = 1 - fit_lsmle$alpha_const)
    fit_lmdpde_r <- LMDPDE_BETA(y, X, Z, gamma_fix = 0, idx_gamma = idx,
                                qoptimal = FALSE, q0 = 1 - fit_lmdpde$alpha_const)
  }

  wald_mle    <- wald_statistic(fit_mle$theta_hat_full, fit_mle$fisher_info, 0, i)
  score_mle   <- score_statistic(fit_mle_r$score, fit_mle_r$fisher_info, i)
  grad_mle    <- gradient_statistic(fit_mle_r$score, fit_mle$theta_hat_full, 0, i)

  wald_lsmle  <- robust_wald_statistic(fit_lsmle$theta_hat, fit_lsmle$Vq$H, fit_lsmle$Vq$J, 0, i)
  score_lsmle <- robust_score_statistic(fit_lsmle_r$score, fit_lsmle_r$Vq$H, fit_lsmle_r$Vq$J, i)
  grad_lsmle  <- robust_gradient_statistic(fit_lsmle_r$score, fit_lsmle$theta_hat, 0,
                                           fit_lsmle_r$Vq$H, fit_lsmle_r$Vq$J, i)

  wald_lmdpde  <- robust_wald_statistic(fit_lmdpde$theta_hat, fit_lmdpde$Vq$H, fit_lmdpde$Vq$J, 0, i)
  score_lmdpde <- robust_score_statistic(fit_lmdpde_r$score, fit_lmdpde_r$Vq$H, fit_lmdpde_r$Vq$J, i)
  grad_lmdpde  <- robust_gradient_statistic(fit_lmdpde_r$score, fit_lmdpde$theta_hat, 0,
                                            fit_lmdpde_r$Vq$H, fit_lmdpde_r$Vq$J, i)

  test_rows[[i]] <- data.frame(
    Parameter = param_names[i],
    Wald_MLE = wald_mle$statistic, Wald_MLE_p = wald_mle$p_value,
    Score_MLE = score_mle$statistic, Score_MLE_p = score_mle$p_value,
    Grad_MLE = grad_mle$statistic, Grad_MLE_p = grad_mle$p_value,
    Wald_LSMLE = wald_lsmle$statistic, Wald_LSMLE_p = wald_lsmle$p_value,
    Score_LSMLE = score_lsmle$statistic, Score_LSMLE_p = score_lsmle$p_value,
    Grad_LSMLE = grad_lsmle$statistic, Grad_LSMLE_p = grad_lsmle$p_value,
    Wald_LMDPDE = wald_lmdpde$statistic, Wald_LMDPDE_p = wald_lmdpde$p_value,
    Score_LMDPDE = score_lmdpde$statistic, Score_LMDPDE_p = score_lmdpde$p_value,
    Grad_LMDPDE = grad_lmdpde$statistic, Grad_LMDPDE_p = grad_lmdpde$p_value
  )
}
test_table <- do.call(rbind, test_rows)
cat("=== Table 4.3 (reduced model): Wald/Rao/Gradient statistics and p-values ===\n")
print(round(test_table[, -1], 4), row.names = FALSE)
cat("\nDone. Compare against Table 4.3 in tese_mestrado/conteudo/05-application.tex.\n")
