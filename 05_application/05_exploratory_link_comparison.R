# SUPPLEMENTARY / EXPLORATORY SCRIPT - does NOT reproduce the dissertation's
# published tables. Scripts 02-04 use GDP on its original scale, matching the
# dissertation (Tables 4.2-4.4); this script instead explores log(GDP) and
# alternative link functions, to document those avenues without conflating
# them with the published results. Kept separate for that reason.
#
# Beta regression application to the HIC dataset, extending aplicacao.R with
# the committee's (banca) suggestions:
#   1. Try alternative link functions for the mean submodel (logit, probit,
#      cloglog, loglog, cauchit) - the rule of thumb raised in the defense is
#      that loglog tends to suit responses concentrated at small values,
#      cloglog suits responses concentrated at large values, and logit/probit
#      (symmetric) suit responses centered in the middle of (0,1). We don't
#      assume which applies here - every link is fit and compared by AIC/BIC,
#      and the comparison itself shows whether the heuristic holds for HIC.
#   2. Add an interaction term (Urb x log(GDP)) to the mean submodel, and
#      formally test whether it improves the fit (likelihood-ratio test)
#      rather than including it by default.
#   3. Report parameter estimates + standard errors, AIC, BIC, and a pseudo-R^2
#      for every candidate, and a LaTeX-ready coefficient table for the final
#      (restricted) model.
#   4. A residual plot with a simulated envelope - the standard diagnostic for
#      beta regression (Ferrari & Cribari-Neto, 2004; Espinheira et al., 2008)
#      - built from this project's own MLE_BETA()/LSMLE_BETA()/LMDPDE_BETA()
#      (no betareg package). This ends up being the most important result
#      below: see the note before Section 5.
#
# Same covariates as aplicacao.R (Urb, GDP for the mean; GDP for the
# precision), but on the log(GDP) scale established in 01_hic_eda.R - unlike
# aplicacao.R, the outlier (row 1, HIC ~ 0.983) is kept in, since the whole
# point of Section 5 is to see how each estimator handles it, not drop it.

rm(list = ls())

if (!dir.exists("../01_functions") || !dir.exists("../02_tests")) {
  stop("Run with the working directory set to 05_application/")
}
invisible(lapply(list.files("../01_functions", pattern = "\\.R$", full.names = TRUE), source))
invisible(lapply(list.files("../02_tests", pattern = "\\.R$", full.names = TRUE), source))

d <- read.table("HICdata.txt", header = TRUE)
y <- d$HIC
n <- length(y)
log_gdp <- log(d$GDP)
outlier_idx <- which.max(y)
dir.create("figures", showWarnings = FALSE)

X_main <- cbind(1, d$Urb, log_gdp)
X_int  <- cbind(X_main, d$Urb * log_gdp)
colnames(X_main) <- c("Intercept", "Urb", "log(GDP)")
colnames(X_int)  <- c(colnames(X_main), "Urb:log(GDP)")
Z <- cbind(1, log_gdp)
colnames(Z) <- c("Intercept", "log(GDP)")

# ---- Helpers: log-likelihood, AIC/BIC, pseudo-R2 ----

# Recomputed directly (not stored by MLE_BETA) - same formula as the internal
# optimized objective in 01_functions/05_mle.R, evaluated at the final estimate
beta_loglik <- function(theta, y, X, Z, linkmu, linkphi = "log") {
  p <- ncol(X); q <- ncol(Z)
  beta <- theta[1:p]; gamma <- theta[(p + 1):(p + q)]
  mu  <- mean_link(beta, y, X, linkmu)$mu
  phi <- precision_link(gamma, Z, linkphi)$phi
  a <- mu * phi; b <- (1 - mu) * phi
  sum(lgamma(phi) - lgamma(a) - lgamma(b) + (a - 1) * log(y) + (b - 1) * log(1 - y))
}

fit_one <- function(y, X, Z, linkmu, label) {
  fit <- tryCatch(MLE_BETA(y, X, Z, linkmu = linkmu), error = function(e) NULL)
  if (is.null(fit)) {
    cat("  [", label, "] did not converge - skipped\n", sep = "")
    return(NULL)
  }
  theta <- fit$theta_hat_full
  ll <- beta_loglik(theta, y, X, Z, linkmu)
  k <- length(theta)
  se <- tryCatch(sqrt(diag(solve(fit$fisher_info))), error = function(e) rep(NA_real_, k))
  mu_hat <- mean_link(theta[1:ncol(X)], y, X, linkmu)$mu
  # NA (not a warning) when X is intercept-only - mu_hat is then a single
  # repeated value and a correlation with it is undefined by construction
  r2_beta <- if (sd(mu_hat) > 0) cor(y, mu_hat)^2 else NA_real_
  list(label = label, linkmu = linkmu, X = X, Z = Z, theta = theta, se = se,
       loglik = ll, aic = -2 * ll + 2 * k, bic = -2 * ll + log(n) * k, k = k,
       r2_beta = r2_beta)
}

# NOTE on pseudo-R2: Cox-Snell/Nagelkerke's formula assumes the log-likelihood
# is <= 0 (true for binary/count GLMs, where likelihood <= 1), so that
# exp((2/n)*loglik_null) <= 1 and the "maximum achievable R2" denominator stays
# positive. The beta distribution is continuous, so its log-likelihood can be
# POSITIVE (density > 1 is normal for a concentrated distribution) - here
# loglik_null ~= 39.7, so exp((2/n)*39.7) > 1 and Nagelkerke's own
# normalization constant comes out NEGATIVE, making the whole statistic
# nonsensical (verified: this is not a code bug, it's a known limitation of
# applying a discrete-outcome pseudo-R2 to a continuous-density model - try
# it and it reproduces for every link/variant below). We report it anyway,
# since it was requested, but flag it, and report the standard beta-regression
# pseudo-R2 instead (squared correlation of y and the fitted mean - the same
# definition the betareg package itself uses) as the one to actually trust.
nagelkerke_r2 <- function(loglik_full, loglik_null, n) {
  r2_cs  <- 1 - exp((2 / n) * (loglik_null - loglik_full))
  max_r2 <- 1 - exp((2 / n) * loglik_null)
  r2_cs / max_r2
}

null_fit <- fit_one(y, matrix(1, n, 1), matrix(1, n, 1), "logit", "null")
ll_null <- null_fit$loglik

# ---- 1. Link function sweep: main effects vs. +interaction, every link ----
links <- c("logit", "probit", "cloglog", "loglog", "cauchit")
results <- list()
for (lk in links) {
  for (variant in c("main", "interaction")) {
    Xv <- if (variant == "main") X_main else X_int
    r <- fit_one(y, Xv, Z, lk, paste(lk, variant, sep = "_"))
    if (!is.null(r)) {
      r$variant <- variant
      r$nagelkerke <- nagelkerke_r2(r$loglik, ll_null, n)
      results[[r$label]] <- r
    }
  }
}

comparison <- do.call(rbind, lapply(results, function(r) data.frame(
  Link = r$linkmu, Interaction = (r$variant == "interaction"),
  LogLik = r$loglik, AIC = r$aic, BIC = r$bic,
  Nagelkerke = r$nagelkerke, R2_beta = r$r2_beta
)))
comparison <- comparison[order(comparison$AIC), ]
cat("\n=== Model comparison across links and interaction (sorted by AIC) ===\n")
cat("(Nagelkerke is included because it was requested, but see the code comment above",
    "nagelkerke_r2() - it is not a valid statistic for a continuous-density model like",
    "beta regression and comes out negative/nonsensical here. R2_beta = cor(y, mu_hat)^2,",
    "the pseudo-R2 the betareg package itself reports, is the one to read.)\n")
print(data.frame(comparison[, 1:2], round(comparison[, 3:7], 4)), row.names = FALSE)

cat("\nHIC's own distribution (from 01_hic_eda.R): mean = 0.225, median = 0.172,",
    "right-skewed - the committee's heuristic (small y -> loglog) predicts",
    "loglog should fit best; check whether the AIC ranking above agrees.\n")

cat("\nCAVEAT: this whole comparison is run with MLE, which Section 5 below shows",
    "is itself badly distorted by the single outlier (row 1). So 'best AIC' here",
    "may just mean 'most flexible at accommodating that one point', not 'genuinely",
    "fits the other 79 countries best' - keep that in mind before over-trusting the",
    "ranking.\n")

# ---- 2. Does adding the interaction help, for the best-fitting link? ----
best_link <- comparison$Link[!comparison$Interaction][1]
r_main <- results[[paste(best_link, "main", sep = "_")]]
r_int  <- results[[paste(best_link, "interaction", sep = "_")]]

lr_stat <- 2 * (r_int$loglik - r_main$loglik)
lr_p    <- 1 - pchisq(lr_stat, df = 1)
cat("\n=== Likelihood-ratio test for the Urb:log(GDP) interaction ===\n")
cat("Best main-effects link:", best_link, "\n")
cat("LR statistic =", round(lr_stat, 4), ", df = 1, p-value =", round(lr_p, 4), "\n")

final <- if (lr_p < 0.05) {
  cat("Interaction is significant at 5% - keeping it in the final model.\n")
  r_int
} else {
  cat("Interaction is NOT significant at 5% - the restricted (main-effects-only)",
      "model is preferred; interaction dropped.\n")
  r_main
}
cat("\nFinal (restricted) model: link =", final$linkmu,
    ", interaction =", (final$variant == "interaction"), "\n")

# ---- 3. Coefficient table for the final (restricted) model ----
coef_names <- c(colnames(final$X), paste0("(phi) ", colnames(Z)))
z_val <- final$theta / final$se
p_val <- 2 * (1 - pnorm(abs(z_val)))
coef_table <- data.frame(
  Parameter = coef_names, Estimate = final$theta, `Std. Error` = final$se,
  `z value` = z_val, `Pr(>|z|)` = p_val, check.names = FALSE
)
cat("\n=== Final (restricted) model - coefficient table ===\n")
print(coef_table, row.names = FALSE, digits = 4)
cat("\nAIC =", round(final$aic, 3), " BIC =", round(final$bic, 3),
    " R2_beta =", round(final$r2_beta, 4), "\n")

# LaTeX-ready version (no external package - just cat() a tabular block)
cat("\n% ---- Copy-paste into LaTeX ----\n")
cat("\\begin{table}[h]\n\\centering\n")
cat("\\caption{Beta regression estimates for the final (restricted) model (link: ",
    final$linkmu, ")}\n", sep = "")
cat("\\begin{tabular}{lrrrr}\n\\hline\n")
cat("Parameter & Estimate & Std. Error & $z$ & $p$-value \\\\\n\\hline\n")
for (i in seq_len(nrow(coef_table))) {
  cat(sprintf("%s & %.4f & %.4f & %.3f & %.4f \\\\\n",
              coef_table$Parameter[i], coef_table$Estimate[i],
              coef_table$`Std. Error`[i], coef_table$`z value`[i],
              coef_table$`Pr(>|z|)`[i]))
}
cat("\\hline\n\\end{tabular}\n\\end{table}\n")

# ---- 4. Residual envelope helper (used in Section 5) ----
# Pearson/ordinary residual: (y - mu_hat) / sqrt(Var(y_hat)), with
# Var(y) = mu(1-mu)/(1+phi) under the beta regression model
# (Ferrari & Cribari-Neto, 2004). Envelope: refit the SAME model to B
# parametric-bootstrap replicates drawn from the fitted model, sort each
# replicate's |residuals|, and take the min/max at each rank as the band -
# the standard construction (Atkinson, 1985; Espinheira et al., 2008).
pearson_resid <- function(y, mu, phi) (y - mu) / sqrt(mu * (1 - mu) / (1 + phi))

residual_envelope <- function(y, X, Z, fit_fun, theta_extract, linkmu, B = 100, seed = 2026) {
  n <- length(y); p <- ncol(X); q <- ncol(Z)
  fit <- fit_fun(y, X, Z)
  theta <- theta_extract(fit)
  mu_hat  <- mean_link(theta[1:p], y, X, linkmu)$mu
  phi_hat <- precision_link(theta[(p + 1):(p + q)], Z, "log")$phi
  obs_resid <- sort(abs(pearson_resid(y, mu_hat, phi_hat)))

  sim_mat <- matrix(NA_real_, n, B)
  set.seed(seed)
  for (b in seq_len(B)) {
    y_sim <- rbeta(n, mu_hat * phi_hat, (1 - mu_hat) * phi_hat)
    fit_b <- tryCatch(fit_fun(y_sim, X, Z), error = function(e) NULL)
    if (is.null(fit_b)) next
    theta_b <- theta_extract(fit_b)
    mu_b  <- mean_link(theta_b[1:p], y_sim, X, linkmu)$mu
    phi_b <- precision_link(theta_b[(p + 1):(p + q)], Z, "log")$phi
    sim_mat[, b] <- sort(abs(pearson_resid(y_sim, mu_b, phi_b)))
  }
  env_lower <- apply(sim_mat, 1, min, na.rm = TRUE)
  env_mid   <- apply(sim_mat, 1, median, na.rm = TRUE)
  env_upper <- apply(sim_mat, 1, max, na.rm = TRUE)
  half_normal_q <- qnorm((seq_len(n) + n - 1 / 8) / (2 * n + 1 / 2))
  n_outside <- sum(obs_resid > env_upper | obs_resid < env_lower)
  list(q = half_normal_q, obs = obs_resid, lower = env_lower, mid = env_mid,
       upper = env_upper, n_outside = n_outside, n = n)
}

plot_envelope <- function(env, title, file) {
  pdf(file, width = 6, height = 5.5)
  plot(env$q, env$obs, xlab = "Half-normal quantiles",
       ylab = "Absolute value of the residuals", pch = 1, col = "black",
       ylim = range(c(env$obs, env$lower, env$upper)), main = title)
  lines(env$q, env$lower, lty = 2, col = "black")
  lines(env$q, env$mid, lty = 1, col = "black")
  lines(env$q, env$upper, lty = 2, col = "black")
  dev.off()
}

# ---- 5. THE KEY RESULT: does the outlier break MLE's residual envelope, and
# do the robust estimators fix it? ----
# Using the simplest, most standard specification (logit, main effects only -
# not the AIC-selected model from Section 2, precisely because that selection
# is itself MLE/outlier-contaminated per the caveat above) so the comparison
# below isolates the ESTIMATOR's robustness, not a particular link's quirks.
cat("\n\n=== Section 5: residual envelope, MLE vs LSMLE vs LMDPDE (all 80 obs, logit) ===\n")

env_mle <- residual_envelope(
  y, X_main, Z, linkmu = "logit",
  fit_fun = function(y, X, Z) MLE_BETA(y, X, Z, linkmu = "logit"),
  theta_extract = function(f) f$theta_hat_full
)
cat("MLE (classical)   - points outside envelope:", env_mle$n_outside, "of", env_mle$n, "\n")

env_lsmle <- residual_envelope(
  y, X_main, Z, linkmu = "logit",
  fit_fun = function(y, X, Z) LSMLE_BETA(y, X, Z, linkmu = "logit"),
  theta_extract = function(f) f$theta_hat
)
cat("LSMLE (robust)    - points outside envelope:", env_lsmle$n_outside, "of", env_lsmle$n, "\n")

env_lmdpde <- residual_envelope(
  y, X_main, Z, linkmu = "logit",
  fit_fun = function(y, X, Z) LMDPDE_BETA(y, X, Z, linkmu = "logit"),
  theta_extract = function(f) f$theta_hat
)
cat("LMDPDE (robust)   - points outside envelope:", env_lmdpde$n_outside, "of", env_lmdpde$n, "\n")

plot_envelope(env_mle, "MLE - residual envelope (n=80, outlier included)",
              "figures/07_residual_envelope_mle.pdf")
plot_envelope(env_lsmle, "LSMLE - residual envelope (n=80, outlier included)",
              "figures/07_residual_envelope_lsmle.pdf")
plot_envelope(env_lmdpde, "LMDPDE - residual envelope (n=80, outlier included)",
              "figures/07_residual_envelope_lmdpde.pdf")

cat("\nIf MLE shows far more violations than LSMLE/LMDPDE here, that's direct, real-data",
    "evidence for this dissertation's central point: a single high-leverage",
    "observation (row 1) distorts the classical fit for the WHOLE sample, while the",
    "robust estimators absorb it without the rest of the fit breaking down.\n")

# ---- 6. Hypothesis tests: classical (MLE) vs robust (LSMLE/LMDPDE) ----
# Same aplicacao.R pattern (per-parameter Wald/Score/Gradient, classical and
# robust), but collected into one table instead of many printed blocks.
# Uses the same baseline design as Section 5 (logit, main effects, all 80
# obs - the well-behaved specification, not the outlier-chasing cauchit +
# interaction one from Sections 1-3) so these tests aren't testing a model
# that's itself a symptom of the contamination problem under study.
cat("\n\n=== Section 6: hypothesis tests per parameter, classical vs robust ===\n")

fit_mle_u    <- MLE_BETA(y, X_main, Z, linkmu = "logit")
fit_lsmle_u  <- LSMLE_BETA(y, X_main, Z, linkmu = "logit")
fit_lmdpde_u <- LMDPDE_BETA(y, X_main, Z, linkmu = "logit")

p_mean <- ncol(X_main); q_prec <- ncol(Z)
param_names <- c(colnames(X_main), paste0("(phi) ", colnames(Z)))
test_rows <- list()

# Only the slopes are tested, not the two intercepts (indices 1 and p_mean+1):
# "intercept = 0" isn't a substantively meaningful hypothesis on the link
# scale (nobody expects mu or phi to collapse to a fixed reference value with
# every covariate at 0), and restricting either intercept to exactly 0 here
# pushes the restricted MLE optimizer into a region where its finite-difference
# gradient blows up (non-finite) for this particular dataset/link combination -
# not worth chasing for hypotheses no one would actually report.
testable_idx <- setdiff(seq_len(p_mean + q_prec), c(1, p_mean + 1))
for (i in testable_idx) {
  is_beta <- i <= p_mean
  idx <- if (is_beta) i else i - p_mean

  if (is_beta) {
    fit_mle_r <- MLE_BETA(y, X_main, Z, linkmu = "logit", beta_fix = 0, idx_beta = idx)
    fit_lsmle_r <- LSMLE_BETA(y, X_main, Z, linkmu = "logit", beta_fix = 0, idx_beta = idx,
                              qoptimal = FALSE, q0 = 1 - fit_lsmle_u$alpha_const)
    fit_lmdpde_r <- LMDPDE_BETA(y, X_main, Z, linkmu = "logit", beta_fix = 0, idx_beta = idx,
                                qoptimal = FALSE, q0 = 1 - fit_lmdpde_u$alpha_const)
  } else {
    fit_mle_r <- MLE_BETA(y, X_main, Z, linkmu = "logit", gamma_fix = 0, idx_gamma = idx)
    fit_lsmle_r <- LSMLE_BETA(y, X_main, Z, linkmu = "logit", gamma_fix = 0, idx_gamma = idx,
                              qoptimal = FALSE, q0 = 1 - fit_lsmle_u$alpha_const)
    fit_lmdpde_r <- LMDPDE_BETA(y, X_main, Z, linkmu = "logit", gamma_fix = 0, idx_gamma = idx,
                                qoptimal = FALSE, q0 = 1 - fit_lmdpde_u$alpha_const)
  }
  test_idx <- i

  wald_mle    <- wald_statistic(fit_mle_u$theta_hat_full, fit_mle_u$fisher_info, 0, test_idx)
  score_mle   <- score_statistic(fit_mle_r$score, fit_mle_r$fisher_info, test_idx)
  grad_mle    <- gradient_statistic(fit_mle_r$score, fit_mle_u$theta_hat_full, 0, test_idx)

  wald_lsmle  <- robust_wald_statistic(fit_lsmle_u$theta_hat, fit_lsmle_u$Vq$H, fit_lsmle_u$Vq$J, 0, test_idx)
  score_lsmle <- robust_score_statistic(fit_lsmle_r$score, fit_lsmle_r$Vq$H, fit_lsmle_r$Vq$J, test_idx)
  grad_lsmle  <- robust_gradient_statistic(fit_lsmle_r$score, fit_lsmle_u$theta_hat, 0,
                                           fit_lsmle_r$Vq$H, fit_lsmle_r$Vq$J, test_idx)

  wald_lmdpde  <- robust_wald_statistic(fit_lmdpde_u$theta_hat, fit_lmdpde_u$Vq$H, fit_lmdpde_u$Vq$J, 0, test_idx)
  score_lmdpde <- robust_score_statistic(fit_lmdpde_r$score, fit_lmdpde_r$Vq$H, fit_lmdpde_r$Vq$J, test_idx)
  grad_lmdpde  <- robust_gradient_statistic(fit_lmdpde_r$score, fit_lmdpde_u$theta_hat, 0,
                                            fit_lmdpde_r$Vq$H, fit_lmdpde_r$Vq$J, test_idx)

  test_rows[[i]] <- data.frame(
    Parameter = param_names[i],
    Wald_MLE = wald_mle$p_value, Score_MLE = score_mle$p_value, Grad_MLE = grad_mle$p_value,
    Wald_LSMLE = wald_lsmle$p_value, Score_LSMLE = score_lsmle$p_value, Grad_LSMLE = grad_lsmle$p_value,
    Wald_LMDPDE = wald_lmdpde$p_value, Score_LMDPDE = score_lmdpde$p_value, Grad_LMDPDE = grad_lmdpde$p_value
  )
}

test_table <- do.call(rbind, test_rows)
cat("\np-values, H0: parameter = 0 (all 80 obs, logit, main effects):\n")
print(cbind(test_table["Parameter"], round(test_table[, -1], 4)), row.names = FALSE)

cat("\nWatch Urb and log(GDP) (the mean submodel's slopes) specifically: if MLE's",
    "p-values disagree sharply with LSMLE/LMDPDE's for those, that's the outlier",
    "changing which covariates look significant depending on which estimator you",
    "trust - directly relevant to what conclusions the dissertation can draw from",
    "this dataset.\n")

cat("\nDone.\n")
