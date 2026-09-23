# Shared engine for the Monte Carlo scenarios (A, B, C).
#
# Run every scenario_*.R script with the working directory set to this folder
# (03_scenarios/), e.g. setwd("path/to/beta-regression-simulations/03_scenarios")
# before sourcing/running a scenario file.

if (!dir.exists("../01_functions") || !dir.exists("../02_tests")) {
  stop(
    "Working directory must be '03_scenarios/' inside the beta-regression-simulations ",
    "repo (found: ", getwd(), "). Use setwd() to that folder first."
  )
}

#**** Only external dependency: MASS::ginv(), used as a fallback when
#     solve() hits a singular matrix. MASS ships with every R installation,
#     so this should never fail; checked once, up front, just in case. ****#
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Missing required package 'MASS'. Install it with install.packages(\"MASS\").")
}

FUNCTION_FILES <- list.files("../01_functions", pattern = "\\.R$", full.names = TRUE, recursive = FALSE)
TEST_FILES     <- list.files("../02_tests", pattern = "\\.R$", full.names = TRUE)
ALL_SOURCE_FILES <- c(FUNCTION_FILES, TEST_FILES)

invisible(lapply(ALL_SOURCE_FILES, source))

#**** Replicates a base design matrix up to n rows ****#
build_X <- function(X_base, n) {
  reps <- n / nrow(X_base)
  if (reps != as.integer(reps)) {
    stop("n = ", n, " must be a multiple of nrow(X_base) = ", nrow(X_base))
  }
  as.matrix(do.call(rbind, replicate(reps, X_base, simplify = FALSE)))
}

#**** Declares one hypothesis restriction. free_idx (which theta positions
#     stay unrestricted) is derived automatically from idx_beta/idx_gamma ****#
make_hypothesis <- function(name, p, q, idx_beta = NULL, beta_fix = NULL,
                            idx_gamma = NULL, gamma_fix = NULL,
                            test_idx, theta0) {
  fixed_idx <- c(idx_beta, if (!is.null(idx_gamma)) p + idx_gamma else NULL)
  list(
    name = name,
    idx_beta = idx_beta, beta_fix = beta_fix,
    idx_gamma = idx_gamma, gamma_fix = gamma_fix,
    test_idx = test_idx, theta0 = theta0,
    free_idx = setdiff(seq_len(p + q), fixed_idx)
  )
}

#**** Fits MLE/LSMLE/LMDPDE under one hypothesis restriction and runs the
#     three classical + three robust tests (Wald, score, gradient) against
#     the corresponding unrestricted fits. Restricted LSMLE/LMDPDE refits
#     always use their OWN unrestricted alpha_const as q0 (see README note
#     on the H1 q0 bug found in the original scripts). ****#
fit_and_test_hypothesis <- function(y, X, Z, hyp, fit_mle, fit_lsmle, fit_lmdpde) {
  fit_mle_h <- try(MLE_BETA(y, X, Z,
    beta_fix = hyp$beta_fix, idx_beta = hyp$idx_beta,
    gamma_fix = hyp$gamma_fix, idx_gamma = hyp$idx_gamma
  ), silent = TRUE)

  fit_lsmle_h <- try(LSMLE_BETA(y, X, Z,
    beta_fix = hyp$beta_fix, idx_beta = hyp$idx_beta,
    gamma_fix = hyp$gamma_fix, idx_gamma = hyp$idx_gamma,
    qoptimal = FALSE, q0 = 1 - fit_lsmle$alpha_const
  ), silent = TRUE)

  fit_lmdpde_h <- try(LMDPDE_BETA(y, X, Z,
    beta_fix = hyp$beta_fix, idx_beta = hyp$idx_beta,
    gamma_fix = hyp$gamma_fix, idx_gamma = hyp$idx_gamma,
    qoptimal = FALSE, q0 = 1 - fit_lmdpde$alpha_const
  ), silent = TRUE)

  fits <- list(mle = fit_mle_h, lsmle = fit_lsmle_h, lmdpde = fit_lmdpde_h)
  if (any(sapply(fits, inherits, "try-error"))) return(list(error = TRUE))

  wald <- list(
    mle    = wald_statistic(fit_mle$theta_hat_full, fit_mle$fisher_info, hyp$theta0, hyp$test_idx),
    lsmle  = robust_wald_statistic(fit_lsmle$theta_hat, fit_lsmle$Vq$H, fit_lsmle$Vq$J, hyp$theta0, hyp$test_idx),
    lmdpde = robust_wald_statistic(fit_lmdpde$theta_hat, fit_lmdpde$Vq$H, fit_lmdpde$Vq$J, hyp$theta0, hyp$test_idx)
  )
  score <- list(
    mle    = score_statistic(fit_mle_h$score, fit_mle_h$fisher_info, hyp$test_idx),
    lsmle  = robust_score_statistic(fit_lsmle_h$score, fit_lsmle_h$Vq$H, fit_lsmle_h$Vq$J, hyp$test_idx),
    lmdpde = robust_score_statistic(fit_lmdpde_h$score, fit_lmdpde_h$Vq$H, fit_lmdpde_h$Vq$J, hyp$test_idx)
  )
  grad <- list(
    mle    = gradient_statistic(fit_mle_h$score, fit_mle$theta_hat_full, hyp$theta0, hyp$test_idx),
    lsmle  = robust_gradient_statistic(fit_lsmle_h$score, fit_lsmle$theta_hat, hyp$theta0, fit_lsmle_h$Vq$H, fit_lsmle_h$Vq$J, hyp$test_idx),
    lmdpde = robust_gradient_statistic(fit_lmdpde_h$score, fit_lmdpde$theta_hat, hyp$theta0, fit_lmdpde_h$Vq$H, fit_lmdpde_h$Vq$J, hyp$test_idx)
  )

  list(error = FALSE, fit = fits, wald = wald, score = score, grad = grad)
}

#**** Assembles one output row (a named list) from a converged replication ****#
build_result_row <- function(i, n, fit_mle, fit_lsmle, fit_lmdpde, hyp_results, hypotheses, param_names) {
  row <- list(
    n = n, sim = i,
    alpha_lsmle = fit_lsmle$alpha_const, alpha_lmdpde = fit_lmdpde$alpha_const
  )

  for (k in seq_along(param_names)) {
    row[[paste0(param_names[k], "_mle")]]    <- fit_mle$theta_hat_full[k]
    row[[paste0(param_names[k], "_lsmle")]]  <- fit_lsmle$theta_hat[k]
    row[[paste0(param_names[k], "_lmdpde")]] <- fit_lmdpde$theta_hat[k]
  }

  for (h_i in seq_along(hypotheses)) {
    hyp <- hypotheses[[h_i]]
    res <- hyp_results[[h_i]]
    tag <- hyp$name

    row[[paste0("alpha_lsmle_", tag)]]  <- res$fit$lsmle$alpha_const
    row[[paste0("alpha_lmdpde_", tag)]] <- res$fit$lmdpde$alpha_const

    for (k in hyp$free_idx) {
      pname <- param_names[k]
      row[[paste0(pname, "_mle_", tag)]]    <- res$fit$mle$theta_hat_full[k]
      row[[paste0(pname, "_lsmle_", tag)]]  <- res$fit$lsmle$theta_hat[k]
      row[[paste0(pname, "_lmdpde_", tag)]] <- res$fit$lmdpde$theta_hat[k]
    }

    for (test_type in c("wald", "score", "grad")) {
      for (est in c("mle", "lsmle", "lmdpde")) {
        stat <- res[[test_type]][[est]]
        row[[paste0("stat_", test_type, "_", tag, "_", est)]] <- stat$statistic
        row[[paste0("pval_", test_type, "_", tag, "_", est)]] <- stat$p_value
      }
    }
  }
  row
}

#**** Extracts a readable message from the first try-error in a named list,
#     e.g. "lsmle: trying to use CRAN without setting a mirror" ****#
first_error_message <- function(named_fits) {
  is_err <- vapply(named_fits, inherits, logical(1), "try-error")
  if (!any(is_err)) return(NULL)
  culprit <- named_fits[is_err][[1]]
  paste0(names(named_fits)[is_err][1], ": ", conditionMessage(attr(culprit, "condition")))
}

#**** Runs one full replication: regenerates data until every fit (unrestricted
#     + every hypothesis) succeeds, then assembles its result row. Gives up
#     after MAX_ATTEMPTS consecutive failures instead of retrying forever, and
#     surfaces the actual underlying error (try() would otherwise swallow it) ****#
MAX_REPLICATION_ATTEMPTS <- 200

run_one_replication <- function(i, cfg) {
  set.seed(cfg$seed_base + i)
  attempt <- 0
  repeat {
    attempt <- attempt + 1
    y <- cfg$generate_data(cfg$X, cfg$Z, cfg$beta, cfg$gamma)

    fit_mle    <- try(MLE_BETA(y, cfg$X, cfg$Z), silent = TRUE)
    fit_lsmle  <- try(LSMLE_BETA(y, cfg$X, cfg$Z), silent = TRUE)
    fit_lmdpde <- try(LMDPDE_BETA(y, cfg$X, cfg$Z), silent = TRUE)
    base_fits <- list(mle = fit_mle, lsmle = fit_lsmle, lmdpde = fit_lmdpde)

    if (any(vapply(base_fits, inherits, logical(1), "try-error"))) {
      if (attempt >= MAX_REPLICATION_ATTEMPTS) {
        stop(
          "Replication ", i, ": ", MAX_REPLICATION_ATTEMPTS,
          " consecutive fit failures, giving up.\nUnderlying error - ",
          first_error_message(base_fits)
        )
      }
      next
    }

    hyp_results <- lapply(cfg$hypotheses, function(h) {
      fit_and_test_hypothesis(y, cfg$X, cfg$Z, h, fit_mle, fit_lsmle, fit_lmdpde)
    })

    if (any(sapply(hyp_results, function(r) isTRUE(r$error)))) {
      if (attempt >= MAX_REPLICATION_ATTEMPTS) {
        stop(
          "Replication ", i, ": ", MAX_REPLICATION_ATTEMPTS,
          " consecutive hypothesis-fit failures, giving up."
        )
      }
      next
    }
    break
  }
  build_result_row(i, cfg$n, fit_mle, fit_lsmle, fit_lmdpde, hyp_results, cfg$hypotheses, cfg$param_names)
}

#**** Runs cfg$N replications (in parallel by default) and returns them as a
#     data.frame - the shared core behind run_scenario() and run_power_study().
#     Progress is only echoed per-replication in sequential mode, since
#     parallel workers can't report back mid-run. ****#
run_replications <- function(cfg) {
  if (isTRUE(cfg$use_parallel)) {
    n_cores <- if (is.null(cfg$n_cores)) max(1, parallel::detectCores() - 1) else cfg$n_cores
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    #**** Each worker re-sources this very file (rather than having the
    #     master hand-pick which objects to clusterExport) so every helper
    #     function/constant defined here - now or later - is automatically
    #     available on workers, with no risk of forgetting one (see the
    #     MAX_REPLICATION_ATTEMPTS bug this replaced). ****#
    wd <- getwd()
    parallel::clusterCall(cl, function(path) { setwd(path); source("00_scenario_utils.R") }, wd)
    parallel::clusterExport(cl, "cfg", envir = environment())

    cat("Running", cfg$N, "replications on", n_cores, "workers...\n")
    results <- parallel::parLapply(cl, seq_len(cfg$N), function(i) run_one_replication(i, cfg))
  } else {
    results <- vector("list", cfg$N)
    for (i in seq_len(cfg$N)) {
      cat("Simulation", i, "of", cfg$N, "| n =", cfg$n, "\n")
      results[[i]] <- run_one_replication(i, cfg)
    }
  }
  do.call(rbind, lapply(results, function(r) as.data.frame(r, check.names = FALSE)))
}

#**** Top-level driver: runs N replications and writes the result table to
#     cfg$output_file ****#
run_scenario <- function(cfg) {
  cat("\n>>> ", cfg$scenario_name, " | n = ", cfg$n, ", N = ", cfg$N, " <<<\n", sep = "")
  final_data <- run_replications(cfg)
  dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
  write.table(final_data, file = cfg$output_file, row.names = FALSE, col.names = TRUE)
  cat("Saved:", cfg$output_file, "\n")
  invisible(final_data)
}

#**** Collapses one grid point's replications into rejection rates: for each
#     pval_* column, the proportion of p-values below alpha (NA-robust - see
#     the na.rm note on nominal_level() in 04_evaluation/01_nominal_level.R).
#     At the null parameter value this rate estimates the test's SIZE; away
#     from the null, it estimates its POWER. ****#
summarize_rejection_rates <- function(rows, param_value, alpha) {
  pval_cols <- grep("^pval_", names(rows), value = TRUE)
  rates <- sapply(rows[pval_cols], function(p) mean(p < alpha, na.rm = TRUE))
  out <- as.data.frame(as.list(rates))
  names(out) <- sub("^pval_", "reject_", pval_cols)
  cbind(param_value = param_value, out)
}

#**** Power-curve driver: sweeps cfg$param_grid (the TRUE value of one
#     parameter used to generate the data), runs cfg$N replications per grid
#     point, and records the rejection rate of cfg$hypothesis at each point.
#     The tested null (cfg$hypothesis$theta0) stays fixed across the whole
#     grid - only the data-generating truth moves. ****#
run_power_study <- function(cfg) {
  cat("\n>>> Power study: ", cfg$study_name, " | n = ", cfg$n, ", N/point = ", cfg$N, " <<<\n", sep = "")

  points <- vector("list", length(cfg$param_grid))
  for (g in seq_along(cfg$param_grid)) {
    value <- cfg$param_grid[g]
    cat("Grid point ", g, "/", length(cfg$param_grid), ": true value = ", value, "\n", sep = "")

    point_cfg <- cfg
    point_cfg[c("beta", "gamma")] <- cfg$param_setter(cfg$beta, cfg$gamma, value)
    point_cfg$hypotheses <- list(cfg$hypothesis) # only the tested hypothesis, nothing else

    rows <- run_replications(point_cfg)
    points[[g]] <- summarize_rejection_rates(rows, value, cfg$alpha)
  }

  power_table <- do.call(rbind, points)
  dir.create(dirname(cfg$output_file), recursive = TRUE, showWarnings = FALSE)
  write.table(power_table, file = cfg$output_file, row.names = FALSE, col.names = TRUE)
  cat("Saved:", cfg$output_file, "\n")
  invisible(power_table)
}
