# beta-regression-simulations

A Monte Carlo simulation study comparing three estimators for beta regression —
maximum likelihood (MLE), the LSMLE (Lq-likelihood / minimum density power
divergence family), and the LMDPDE — under three classical hypothesis-testing
procedures: Wald, Score (Rao), and Gradient. The main question we're after is
practical: when the data are clean, how do these tests behave, and when a
small fraction of the sample is contaminated, which estimator/test
combinations still hold up?

Beta regression follows the Ferrari–Cribari-Neto parametrization: a mean
submodel for μ and a separate precision submodel for φ, each with its own
link function. MLE is the standard choice but is well known to be sensitive
to outliers; LSMLE and LMDPDE are robust alternatives built on density power
divergences, with a tuning parameter that trades off efficiency against
robustness.

## What's in here

- **`01_functions/`** — the estimators themselves: link functions, initial
  values, the MLE score/Fisher information, and the LSMLE/LMDPDE fitting
  routines.
- **`02_tests/`** — the six test statistics: classical and robust versions of
  Wald, Score, and Gradient.
- **`03_scenarios/`** — the simulation engine (`00_scenario_utils.R`) plus
  three data-generating scenarios (A, B, C), each run with and without
  contamination. Contamination replaces a small fraction of the most extreme
  observations with draws pulled toward a fixed target, simulating the kind
  of localized data corruption a robust estimator is supposed to withstand.
  This folder also has the power-curve studies, which sweep the true
  parameter value across a grid and track each test's rejection rate.
- **`04_evaluation/`** — scripts and write-ups that turn the raw simulation
  output into something readable: nominal-level tables, quantile checks, and
  power-curve plots. It also includes a deeper dive into a specific finding —
  under contamination, the Score test's power can *drop* as the true
  parameter moves further from the null, instead of rising monotonically.
  `score_dip_derivation.pdf` works out why, from first principles.
- **`05_application/`** — a real-data example applying these models outside
  of simulation.

## A finding worth flagging

One of the more interesting results from this project: under contamination,
the classical Score test's power curve is non-monotonic — it rises, peaks,
and then falls as the true parameter moves further from the null. This
happens because a couple of high-leverage contaminated points pull the score
in the opposite direction from the rest of the sample, and past a certain
distance from the null the two effects partially cancel out. The robust
LSMLE/LMDPDE-based Score test shows the same dip at small sample sizes, but
it fades away as n grows, since the density-power weighting increasingly
discounts the contaminated points as they diverge from what the model
predicts. The classical MLE-based test doesn't have that luxury — its size
is already broken under contamination regardless of n. The full derivation,
including why the Gradient test doesn't share this weakness, is in
`04_evaluation/score_dip_derivation.pdf`.

## Running it

Everything is base R (plus `MASS` for generalized inverses and `parallel`
for the Monte Carlo loops) — no `ggplot2`, no modeling packages beyond what's
implemented here. Scripts assume the working directory is set to the folder
they live in (e.g. `03_scenarios/` before running a scenario script), since
they use relative paths to `source()` the functions they depend on.
