# Pelvic Curvature Location Analysis (`curvatureLocation.R`)

This script asks whether developmental environment (Aquatic vs. Terrestrial) predicts where along the body peak dorsal curvature tends to occur during a locomotor cycle — Anterior vs. Posterior. It's a binomial GLM at the clip level: environment as the only fixed effect, no random effect, followed by a predicted-probability plot and a power analysis.

## What it does, end to end

1. Loads long-form (cycle-level) curvature-location data and collapses `curv_loc_bin` into Anterior (bins 2–3) or Posterior (bins 4–5); any other bin value is dropped as a sparse edge bin.
2. Aggregates to one row per clip: `posterior_cycles`, `anterior_cycles`, `total_cycles`, and `posterior_prop`.
3. Deduplicates to exactly one clip per fish (see "One clip per fish" below).
4. Fits `cbind(posterior_cycles, anterior_cycles) ~ environment` as a binomial GLM.
5. Reports the model summary, a Type II analysis of deviance (`car::Anova`), and an odds-ratio table with 95% CIs.
6. Computes and plots the model's predicted probability of posterior peak curvature by environment, with 95% CI error bars, and saves the plot.
7. Runs base R's diagnostic plots for the model.
8. Runs a parametric-bootstrap power analysis for the environment effect.

## One clip per fish

Every fish contributes exactly one clip —  This makes "clip" and "fish" the same unit of replication: the response variable (posterior vs. anterior cycles per clip) treats the clip as the unit of replication, and because no fish contributed more than one clip, there's no pseudoreplication to correct for at the fish level. It also means a random effect of fish identity isn't included — with exactly one observation per fish there are no repeated measures within individuals for a `(1 | fishID)` term to account for, so it wouldn't be identifiable even if added.

## Input

| File | Description |
|---|---|
| Long-form curvature-location data | One row per cycle, with `curv_loc_bin`, `clip`, `fishID`, `environment`, `training`. `environment`/`training` are stored as text labels ("Aquatic"/"Terrestrial", "Trained"/"Untrained"), so they're converted straight to factors without relabeling from numeric codes. |

Input path and `out_dir` are set near the top of the script — update them for your file locations.

## Statistical design

```
cbind(posterior_cycles, anterior_cycles) ~ environment
```

A binomial GLM with known trial counts per clip (`total_cycles`), one row per fish.

**Predicted-probability CI**: `predict(..., type = "response", se.fit = TRUE)` gives the fitted probability and a standard error on the probability scale. Since a probability-scale SE isn't valid to use directly inside a logit-scale interval, it's converted to the logit scale first via the standard delta-method transform (dividing by p(1−p), the derivative of the logit link) before building the CI as `plogis(qlogis(p) ± 1.96 × SE_link)` — this keeps the interval properly bounded to [0, 1] and asymmetric near the edges, and is equivalent to what `predict(type = "link", se.fit = TRUE)` would give directly.

## Result

n = 24 (one row per fish).

| Term | Estimate | SE | z | p |
|---|---|---|---|---|
| (Intercept) | -0.571 | 0.347 | -1.644 | 0.100 |
| environmentTerrestrial | 1.041 | 0.478 | 2.176 | **0.030** |

Type II analysis of deviance: environment, LR χ² = 4.897, df = 1, **p = 0.027**.

Odds ratio (Terrestrial vs. Aquatic): **2.83 [1.11, 7.23]** (Wald-based CI, matching the odds-ratio table's convention; a profile-likelihood CI gives a nearly identical **[1.13, 7.40]**).

Predicted probability of posterior peak curvature: Aquatic **36.1% [22.3%, 52.7%]**, Terrestrial **61.5% [45.6%, 75.3%]**.

## Power analysis

A frequentist parametric bootstrap, testing detectability at the model's own fitted effect size (not swept across a range). For each of 1,000 simulations: draw new `posterior_cycles` counts per fish from `rbinom(n, size = total_cycles, prob = fitted_probability)`, holding each fish's real `total_cycles` fixed, refit the same GLM on the simulated counts, and record whether the environment term's p-value clears .05. This answers "if the true effect is what this fitted model says — an odds ratio of ~2.83 — how often would a fresh sample of this same size (n=24) detect it?"

Result: **60.5% power [57.4%, 63.5%]**. Roughly 3 in 5 similarly-sized studies would find a significant environment effect at this sample size and effect size — worth keeping in mind alongside the p = 0.030 result, especially if you're weighing whether to collect more fish before treating a future non-significant result as a true null.

## Output

Written to `out_dir`:

| File | Contents |
|---|---|
| `predictedProb_environment.svg` | Predicted-probability bar plot with 95% CI error bars |
| `curvatureLocation_environmentModel_oddsRatios.csv` | Odds-ratio table (estimate, SE, z, p, odds ratio, 95% CI) |
| `curvatureLocation_environmentModel_predictedProbabilities.csv` | Predicted probability by environment, with SE and 95% CI |
| `curvatureLocationPower_deduplicated.csv` | Power summary for the environment effect (one row) |

## Dependencies

R packages: `readr`, `dplyr`, `car`, `ggplot2`.

Install any missing packages with:
```r
install.packages(c("readr", "dplyr", "car", "ggplot2"))
```

## Running it

No known blocking bugs. Update the input path and `out_dir` near the top if your file locations differ, then run top-to-bottom. The model fit and plot are effectively instant; the power analysis (1,000 `glm()` refits on 24 rows) adds a couple of seconds at most.
