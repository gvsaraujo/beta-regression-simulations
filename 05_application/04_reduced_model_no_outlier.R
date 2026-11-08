# Reduced model (Urb excluded from the precision submodel), MLE only, with
# Observation #1 (the influential outlier, HIC ~ 0.983) removed from the
# dataset. Reproduces Table 4.4 of the dissertation: with the outlier gone,
# the classical MLE's conclusions should move closer to the robust
# estimators' conclusions from Table 4.3 (fit on the full, contaminated
# data) - direct evidence that the robust estimators were already
# compensating for this single observation.

if (!dir.exists("../01_functions") || !dir.exists("../02_tests")) {
  stop("Run with the working directory set to 05_application/")
}
invisible(lapply(list.files("../01_functions", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("../02_tests", pattern = "\\.R$", full.names = TRUE), source))

d <- read.table("HICdata.txt", header = TRUE)[-1, ] # drop Observation #1
y <- d$HIC
X <- cbind(1, d$Urb, d$GDP)
Z <- cbind(1, d$GDP)
param_names <- c("beta1 (Int.)", "beta2 (Urb)", "beta3 (GDP)",
                 "gamma1 (Int.)", "gamma2 (GDP)")

fit_mle <- MLE_BETA(y, X, Z)
se_mle <- sqrt(diag(solve(fit_mle$fisher_info)))
cat("=== Table 4.4 (reduced model, Observation #1 removed): MLE estimates (SE) ===\n")
print(data.frame(Parameter = param_names, Estimate = round(fit_mle$theta_hat_full, 4),
                 SE = round(se_mle, 4)), row.names = FALSE)

test_rows <- list()
for (i in seq_along(param_names)) {
  is_beta <- i <= 3
  idx <- if (is_beta) i else i - 3
  fit_mle_r <- if (is_beta) MLE_BETA(y, X, Z, beta_fix = 0, idx_beta = idx) else
    MLE_BETA(y, X, Z, gamma_fix = 0, idx_gamma = idx)

  wald  <- wald_statistic(fit_mle$theta_hat_full, fit_mle$fisher_info, 0, i)
  score <- score_statistic(fit_mle_r$score, fit_mle_r$fisher_info, i)
  grad  <- gradient_statistic(fit_mle_r$score, fit_mle$theta_hat_full, 0, i)

  test_rows[[i]] <- data.frame(
    Parameter = param_names[i],
    Wald = wald$statistic, Wald_p = wald$p_value,
    Rao = score$statistic, Rao_p = score$p_value,
    Gradient = grad$statistic, Gradient_p = grad$p_value
  )
}
test_table <- do.call(rbind, test_rows)
cat("\n=== Table 4.4: Wald/Rao/Gradient statistics and p-values (MLE) ===\n")
print(round(test_table[, -1], 4), row.names = FALSE)
cat("\nDone. Compare against Table 4.4 in tese_mestrado/conteudo/05-application.tex.\n")
