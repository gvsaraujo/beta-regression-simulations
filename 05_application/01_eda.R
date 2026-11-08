# Exploratory analysis of the HIC dataset used in the real-data application
# (Queiroz-Maluf beta regression example). 80 countries, 3 variables:
#   Urb - urbanization rate (proportion, already in (0,1))
#   GDP - GDP-related measure (raw scale, heavily right-skewed)
#   HIC - the response: bounded in (0,1), the natural target for a beta
#         regression, since it can't be modeled well with an unbounded
#         linear model (a fitted line would eventually predict values
#         outside (0,1) for extreme covariate values)
# Base R graphics only, consistent with the rest of the project - no ggplot2.

rm(list = ls())

d <- read.table("HICdata.txt", header = TRUE)
dir.create("figures", showWarnings = FALSE)

# The one clear outlier in this dataset: HIC ~ 0.983, right up against the
# (0,1) boundary beta regression requires - flagged once here and highlighted
# ONLY in plots involving HIC (never in plots that are just about Urb/GDP).
outlier_idx <- which.max(d$HIC)
cat("\nOutlier (max HIC):\n")
print(d[outlier_idx, ])

# No color - regular points are open circles, bars are faint transparent
# black, and the outlier is always a solid filled black dot
fill_col <- rgb(0, 0, 0, alpha = 0.08)

# ---- 1. Basic structure ----
cat("Dimensions:", nrow(d), "observations,", ncol(d), "variables\n")
str(d)
cat("\nMissing values per column:\n")
print(colSums(is.na(d)))

# ---- 2. Summary statistics ----
cat("\nSummary statistics:\n")
print(summary(d))

# Mean, median, skewness, kurtosis for every variable (committee request),
# log(GDP) included since it's the scale used everywhere from here on down
log_gdp <- log(d$GDP)
skewness <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurtosis <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3 # excess kurtosis: 0 = normal-like tails

desc_stats <- data.frame(
  Variable = c("Urb", "GDP", "log(GDP)", "HIC"),
  Mean     = c(mean(d$Urb), mean(d$GDP), mean(log_gdp), mean(d$HIC)),
  Median   = c(median(d$Urb), median(d$GDP), median(log_gdp), median(d$HIC)),
  SD       = c(sd(d$Urb), sd(d$GDP), sd(log_gdp), sd(d$HIC)),
  Skewness = c(skewness(d$Urb), skewness(d$GDP), skewness(log_gdp), skewness(d$HIC)),
  Kurtosis = c(kurtosis(d$Urb), kurtosis(d$GDP), kurtosis(log_gdp), kurtosis(d$HIC))
)
cat("\nDescriptive statistics (mean, median, SD, skewness, excess kurtosis):\n")
print(round(desc_stats[, -1], 4), row.names = desc_stats$Variable)

# ---- 3. GDP is heavily right-skewed - work in log(GDP) from here on ----
# skewness(GDP) ~ 3 (one country at ~200 against a median of ~25);
# log(GDP) is much closer to symmetric. EVERY later step (correlation,
# distributions, scatterplots, dispersion grouping, pairs) uses log_gdp
# instead of raw GDP - this is the only place the raw scale is shown, as the
# before/after justification for that choice.
pdf("figures/02_log_gdp.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
hist(d$GDP, main = "", xlab = "GDP", col = fill_col, border = "black")
hist(log_gdp, main = "", xlab = "log(GDP)", col = fill_col, border = "black")
dev.off()

# ---- 4. HIC must sit strictly inside (0,1) for beta regression to apply ----
cat("\nHIC range:", paste(round(range(d$HIC), 4), collapse = " to "), "\n")
cat("Any HIC <= 0 or >= 1?", any(d$HIC <= 0 | d$HIC >= 1), "\n")

# ---- 5. Correlation matrix (Urb, log(GDP), HIC) ----
d_model <- data.frame(Urb = d$Urb, `log(GDP)` = log_gdp, HIC = d$HIC, check.names = FALSE)
cat("\nCorrelation matrix:\n")
print(round(cor(d_model), 3))
# Urb and log(GDP) are only weakly correlated with each other (~0.13 - note
# this is higher than raw GDP's ~-0.03, since log() reshapes the relationship
# with Urb too, not just GDP's own distribution), so multicollinearity
# between them isn't a concern; both correlate moderately with HIC in the
# same direction.

# ---- 6. Distributions: histogram + boxplot for each variable ----
pdf("figures/01_distributions.pdf", width = 10, height = 7)
par(mfrow = c(2, 3))
for (v in names(d_model)) {
  hist(d_model[[v]], main = "", xlab = v, col = fill_col, border = "black")
  if (v == "HIC") abline(v = d_model[[v]][outlier_idx], col = "black", lwd = 2, lty = 2)
}
for (v in names(d_model)) {
  boxplot(d_model[[v]], main = "", ylab = v, col = fill_col)
  if (v == "HIC") points(1, d_model[[v]][outlier_idx], col = "black", pch = 19, cex = 1.0)
}
dev.off()

# ---- 7. Response vs covariates ----
pdf("figures/03_scatterplots.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
plot(d$Urb, d$HIC, main = "", xlab = "Urb", ylab = "HIC", pch = 1, col = "black")
lines(lowess(d$Urb, d$HIC), col = "black", lwd = 1.5)
points(d$Urb[outlier_idx], d$HIC[outlier_idx], col = "black", pch = 19, cex = 1.2)
plot(log_gdp, d$HIC, main = "", xlab = "log(GDP)", ylab = "HIC", pch = 1, col = "black")
lines(lowess(log_gdp, d$HIC), col = "black", lwd = 1.5)
points(log_gdp[outlier_idx], d$HIC[outlier_idx], col = "black", pch = 19, cex = 1.2)
dev.off()

# ---- 8. Response on the logit scale ----
# Beta regression (logit link) models logit(mu) as a linear combination of
# the betas - so it's the logit(HIC)-vs-covariate relationship that should
# look roughly linear, not the raw HIC-vs-covariate one.
logit <- function(p) log(p / (1 - p))
pdf("figures/04_logit_scale.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
plot(d$Urb, logit(d$HIC), main = "", xlab = "Urb", ylab = "logit(HIC)", pch = 1, col = "black")
lines(lowess(d$Urb, logit(d$HIC)), col = "black", lwd = 1.5)
points(d$Urb[outlier_idx], logit(d$HIC[outlier_idx]), col = "black", pch = 19, cex = 1.2)
plot(log_gdp, logit(d$HIC), main = "", xlab = "log(GDP)", ylab = "logit(HIC)", pch = 1, col = "black")
lines(lowess(log_gdp, logit(d$HIC)), col = "black", lwd = 1.5)
points(log_gdp[outlier_idx], logit(d$HIC[outlier_idx]), col = "black", pch = 19, cex = 1.2)
dev.off()

# ---- 9. Dispersion pattern: does HIC's spread change with the covariates? ----
# Beta regression fits a separate precision submodel (phi), on top of the
# mean submodel - it's worth checking informally whether HIC looks more/less
# variable across the covariate range, which is exactly what phi would pick
# up. Split log(GDP) into terciles and compare HIC's spread across groups
# (same grouping as raw-GDP terciles would give, since log() is monotonic -
# using log_gdp here just keeps the scale consistent with everything else).
gdp_tercile <- cut(log_gdp, quantile(log_gdp, probs = c(0, 1/3, 2/3, 1)),
                   include.lowest = TRUE, labels = c("low", "mid", "high"))
cat("\nHIC spread (sd) by log(GDP) tercile:\n")
print(tapply(d$HIC, gdp_tercile, sd))

pdf("figures/05_dispersion_by_gdp_tercile.pdf", width = 5, height = 5)
boxplot(d$HIC ~ gdp_tercile, main = "", xlab = "log(GDP) tercile", ylab = "HIC", col = fill_col)
points(as.integer(gdp_tercile[outlier_idx]), d$HIC[outlier_idx], col = "black", pch = 19, cex = 1.0)
dev.off()

# ---- 10. Pairwise scatterplot matrix (Urb, log(GDP), HIC) ----
# Custom panel: highlights the outlier only on panels where one of the two
# axes IS HIC (so the Urb-vs-log(GDP) panel is left plain, untouched)
panel_hic_highlight <- function(x, y, ...) {
  points(x, y, pch = 1, col = "black")
  if (identical(x, d_model$HIC) || identical(y, d_model$HIC)) {
    points(x[outlier_idx], y[outlier_idx], pch = 19, col = "black", cex = 1.2)
  }
}
pdf("figures/06_pairs.pdf", width = 7, height = 7)
pairs(d_model, main = "", panel = panel_hic_highlight)
dev.off()

# ---- 11. Two figures curated for the dissertation (tese_mestrado) ----
# A minimal pair chosen out of the figures above, to avoid overloading the
# Application chapter with every exploratory plot: histograms justifying a
# log(GDP) view of the skewed covariate, and the HIC vs log(GDP) scatter
# with the influential observation highlighted (ties to Figure 2.1).
dir.create("../tese_mestrado/figures", showWarnings = FALSE, recursive = TRUE)

pdf("../tese_mestrado/figures/Application_HIC_GDP_Distributions.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
hist(d$GDP, main = "", xlab = "GDP", col = fill_col, border = "black")
hist(log_gdp, main = "", xlab = "log(GDP)", col = fill_col, border = "black")
dev.off()

pdf("../tese_mestrado/figures/Application_HIC_vs_logGDP.pdf", width = 5.5, height = 4.5)
plot(log_gdp, d$HIC, xlab = "log(GDP)", ylab = "HIC", pch = 1, col = "black")
lines(lowess(log_gdp, d$HIC), col = "black", lwd = 1.5)
points(log_gdp[outlier_idx], d$HIC[outlier_idx], col = "black", pch = 19, cex = 1.3)
dev.off()

cat("\nAll figures saved under 05_application/figures/ and in tese_mestrado/figures/.\n")
