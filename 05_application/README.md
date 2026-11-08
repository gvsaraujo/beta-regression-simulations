# HIC application

Real-data application of the robust beta regression tests (Chapter 4 of the
dissertation) to the health insurance coverage (HIC) dataset.

## Data

`HICdata.txt`: 80 municipalities in the state of São Paulo, Brazil, 2010.
Three columns: `Urb` (percentage of the population in the urban zone), `GDP`
(per capita gross domestic product), `HIC` (health insurance coverage
index, the response). Originally collected by the Institute of Applied
Economic Research (IPEA) and used by Maluf, Ferrari & Queiroz (2025),
*"Robust beta regression through the logit transformation"*, Metrika, 88(1),
61-81. Data and original code: <https://github.com/ffqueiroz/RobustBetareg>.

Observation #1 (HIC ≈ 0.983) is an influential outlier discussed throughout
this chapter.

## How to reproduce

Run from this directory (`05_application/`), in order; each script is
self-contained (sources `../01_functions` and `../02_tests`, reads
`HICdata.txt`) and can also be run on its own.

| Script | Reproduces |
|---|---|
| `01_eda.R` | Table 4.1 (descriptive statistics) and the two application figures used in the dissertation (saved directly to `../tese_mestrado/figures/`), plus additional exploratory plots in `figures/` |
| `02_unrestricted_model.R` | Table 4.2 (unrestricted-model estimates and per-parameter Wald/Rao/Gradient tests, all 80 observations) |
| `03_reduced_model.R` | Table 4.3 (reduced model, Urb dropped from the precision submodel, all 80 observations) |
| `04_reduced_model_no_outlier.R` | Table 4.4 (reduced model, MLE only, Observation #1 removed) |
| `05_exploratory_link_comparison.R` | Supplementary only - **does not** reproduce the published tables. Explores log(GDP) and alternative link functions instead of the GDP scale used in `02`-`04`; kept separate for that reason. |

All four main scripts (`02`-`04`) use GDP on its original scale, matching
`\citet{maluf2025robust}` and the dissertation's model equations.

## Requirements

R, with `MASS` (used as a fallback by the project's own estimator
functions). No other external packages.

## Notes

- `old/` holds earlier, superseded drafts of this analysis (kept for
  history, not part of the reproduction path above).
- The LSMLE and LMDPDE columns of Tables 4.2 and 4.3 were corrected during
  review: an earlier version of the dissertation had the LMDPDE column
  accidentally duplicated from LSMLE. The numbers produced by `02_` and
  `03_` here are the verified, correct ones.
