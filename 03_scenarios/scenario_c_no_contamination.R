# Scenario C: unbounded beta densities, varying precision. No contamination.
# beta = (-3, 7.5), gamma = (1, 2); mu in (0.05, 0.98), phi in (2.7, 20.1).
#
# Unlike A/B, the precision submodel here has its own covariate (gamma2). It
# reuses the same design matrix as the mean submodel, so Z = X.

rm(list = ls())

# Run with the working directory set to 03_scenarios/
source("00_scenario_utils.R")

X_matrix <- read.table("data/X_matrix_40.txt", header = TRUE)

sample_sizes <- c(40, 80, 160, 320, 1200)
N <- 10000
beta  <- c(-3, 7.5)
gamma <- c(1, 2)
p <- length(beta); q <- length(gamma)
param_names <- c("beta1", "beta2", "gamma1", "gamma2")

hypotheses <- list(
  make_hypothesis("h1", p, q, idx_beta = 2, beta_fix = beta[2],
                   test_idx = 2, theta0 = beta[2]),
  make_hypothesis("h2", p, q, idx_gamma = 1, gamma_fix = gamma[1],
                   test_idx = 3, theta0 = gamma[1]),
  make_hypothesis("h3", p, q, idx_beta = c(1, 2), beta_fix = beta,
                   test_idx = c(1, 2), theta0 = beta),
  make_hypothesis("h4", p, q, idx_beta = 2, beta_fix = beta[2],
                   idx_gamma = 2, gamma_fix = gamma[2],
                   test_idx = c(2, 4), theta0 = c(beta[2], gamma[2])),
  make_hypothesis("h5", p, q, idx_beta = c(1, 2), beta_fix = beta,
                   idx_gamma = 2, gamma_fix = gamma[2],
                   test_idx = c(1, 2, 4), theta0 = c(beta, gamma[2])),
  make_hypothesis("h6", p, q, idx_gamma = 2, gamma_fix = gamma[2],
                   test_idx = 4, theta0 = gamma[2])
)

generate_data <- function(X, Z, beta, gamma) {
  mu <- exp(X %*% beta) / (1 + exp(X %*% beta))
  phi <- exp(Z %*% gamma)
  rbeta(nrow(X), mu * phi, (1 - mu) * phi)
}

for (n in sample_sizes) {
  X <- build_X(X_matrix, n)
  Z <- X # precision submodel reuses the mean submodel's covariates

  cfg <- list(
    scenario_name = "Scenario C - no contamination",
    n = n, N = N, beta = beta, gamma = gamma,
    X = X, Z = Z,
    param_names = param_names, hypotheses = hypotheses,
    generate_data = generate_data,
    seed_base = 123,
    use_parallel = TRUE, n_cores = NULL,
    output_file = paste0("Results/scenario_c_no_contamination_n", n, ".txt")
  )

  run_scenario(cfg)
}
