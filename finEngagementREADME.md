# Fin Engagement Analysis Pipeline

This R script (`finEngagementAnnotated.R`) asks whether the fins were engaged for pooled with in each environment at all during a trial, not which kind of fin use occurred, just whether any occurred. 
It fits a frequentist binomial GLMM (`lme4::glmer`), runs a fairly thorough diagnostic suite, and finishes with a simulation-based power analysis. 

## What it does, end to end

1. Loads the treadmill dataset and collapses the 3-level `finUse` column into a binary `finEngaged` (0 = "None", 1 = "Asymmetrical" or "Symmetrical").
2. Fits `finEngaged ~ environment * flowSpeed + (1 | fishID)` as a binomial GLMM.
3. Saves the model summary to a text file, and a tidy fixed-effects table to CSV (twice — once plain, once with odds ratios added).
4. Computes model-predicted probabilities of fin engagement across environment × flowSpeed, with 95% confidence intervals built correctly on the link scale, and saves them to CSV.
5. Computes observed (raw, non-model) proportions of fin engagement per environment × flowSpeed cell, and plots them alongside the model's fitted curve and CI ribbon, saving the plot as an SVG.
6. Runs a GLMM diagnostic suite: DHARMa simulated-residual tests (uniformity, dispersion, outliers), `performance` package checks (singularity, collinearity, an overall model-check panel), and per-fish Cook's distance via `influence.ME`.
7. Runs a simulation-based power analysis for the environment, flowSpeed, and environment:flowSpeed effects.

## A note on model scope

The fixed-effect formula here is just `environment * flowSpeed` — main effects of environment and flowSpeed, plus their interaction. 
This model purposefully pools data within environments together. 

## Input

| File | Description |
|---|---|
| `combinedWalkDataSend.csv` | One row per trial/timepoint, repeated per `fishID`, with a `finUse` column. |

Default input path (edit near the top of the script to change):
```
combinedWalkDataSend.csv
```

## Statistical design

```
finEngaged ~ environment * flowSpeed + (1 | fishID)
```

`finEngaged` starts as `ifelse(finUse == "None", 0, 1)`, then gets converted to a 2-level factor (`"NoFinUse"`/`"FinEngaged"`) for readability, then converted back to numeric (`as.numeric(finEngaged) - 1`) right in the model formula — a roundabout path, but since `factor()` assigns integer codes in the order of its `levels` argument, this correctly recovers the original 0/1 coding.

**Predicted probabilities** are computed on the link (logit) scale first (`type = "link"`, population-level via `re.form = NA`), then back-transformed through `plogis()` for both the point estimate and the 95% CI bounds — this keeps the CI properly bounded to [0, 1] and asymmetric near the edges, which is the statistically preferable approach over building a CI directly on the probability scale. A delta-method approximate SE on the probability scale (`se_prob`) is also computed and saved, but isn't actually used for the reported CI — it's informational only in the output CSV.

**Diagnostics**: DHARMa's `simulateResiduals()` approach simulates many replicate datasets from the fitted model and checks where each real observation falls within its simulated distribution — a way of getting interpretable residual diagnostics for a binomial GLMM, where raw/Pearson residuals aren't straightforward to read the way OLS residuals are. `performance::check_model()` bundles several more checks (including collinearity and singularity) into one panel. `influence.ME` refits the model leaving each fish out in turn and computes Cook's distance per fish, flagging any single fish that disproportionately drives the fitted coefficients (against the common 4/n rule-of-thumb threshold, shown as the plot's red dashed line).

**Important**: none of the diagnostic output (the DHARMa plots, the uniformity/dispersion/outlier test results, `check_model()`'s panel, the Cook's-distance plot) is written to a file anywhere — it's all printed or plotted to the active R graphics device / console only. If you want to keep any of it, you'd need to wrap those calls in a graphics-device call (`pdf()`/`svg()` … `dev.off()`) or capture text output the way the model summary already is (via `sink()`).

## Power analysis

This is a *retrospective/self-consistency* power check — "if the world really looks like what this fitted model says, how often would each effect be detected in a fresh sample of the same size?" — tested only at the model's own observed effect size, not swept across multiple hypothetical sizes.

The mechanism is a standard frequentist parametric bootstrap, using `lme4`'s own `simulate()` method for `glmer` models. It fixes the fixed effects at their fitted (MLE) point estimates and draws fresh random-effect values (from the fitted variance components) plus fresh binomial outcomes for each simulated dataset. That means it doesn't propagate parameter *uncertainty* into the power estimate — only random-effect and residual/binomial sampling variability around the fitted point estimate.

Mechanically:
1. Generate all 1,000 simulated response vectors up front in one call (`simulate(mod_fin, nsim = 1000, seed = 123)`).
2. For each simulated replicate, refit the same GLMM (with the `bobyqa` optimizer explicitly set), and skip it — excluding it from the power estimate's denominator — if the refit either fails outright or raises an lme4 convergence warning.
3. For surviving fits, record whether each of the environment, flowSpeed, and environment:flowSpeed terms' p-value cleared p < .05.
4. Checkpoint the running results to CSV every 100 simulations.

The final table reports, per effect, the proportion of valid simulations where it was "detected," with a binomial exact confidence interval.

## Output

Written to the same directory as the input file:

| File | Contents |
|---|---|
| `finEngagement_model_summary.txt` | Full `summary(mod_fin)` output, captured via `sink()` |
| `finEngagementEffects.csv` | Tidy fixed-effects table (estimate, SE, statistic, p-value) |
| `finEngagementEffects_OR.csv` | Same table, plus an `odds_ratio` column — a superset of the file above, not a separate analysis |
| `finEngagementProbabilities.csv` | Model-predicted probability of fin engagement per environment × flowSpeed, with link-scale-derived 95% CI and the (unused) delta-method `se_prob` column |
| `FinEngagement_ModelObserved.svg` | Plot: model fit (line + CI ribbon) vs. observed proportions (points), by environment across flow speed |
| `finEngagementPower_partial.csv` | Checkpoint file, updated every 100 power-simulation iterations |
| `finEngagementPower.csv` | Final power summary — one row per effect (environment, flowSpeed, environment:flowSpeed), with the estimated detection probability |

Default output directory (edit near the top of the script to change):
```
dataOutputPolyp
```

## Dependencies

R packages: `dplyr`, `readr`, `lme4`, `broom.mixed`, `ggplot2`, `DHARMa`, `performance`, `car`, `influence.ME`.

`car` is loaded in the diagnostics section but no `car::` function is called directly anywhere in the script — `performance::check_collinearity()` may use it internally as a dependency, but it isn't invoked explicitly here.

Install any missing packages with:
```r
install.packages(c("dplyr", "readr", "lme4", "broom.mixed", "ggplot2",
                    "DHARMa", "performance", "car", "influence.ME"))
```

## Running it

No known blocking bugs, so it's safe to `source()` or `Rscript` top-to-bottom.

1. Update the input path near the top and the `out_dir` path if your file locations differ.
2. The main GLMM fit and diagnostics should be quick — this is a relatively lightweight model.
3. The power-analysis loop (plain `glmer` refits, not Stan/MCMC) is also fast, so the whole script should complete in a reasonable amount of time.

## Notes / things to double-check when reading results

- Confirm the model's scope (no `training` term) matches what you actually want tested — see "A note on model scope" above.
- The diagnostic plots and test results aren't saved anywhere; if you're running this non-interactively (e.g. via `Rscript`), you won't see them at all unless you add a graphics-device wrapper.
- `se_prob` in `finEngagementProbabilities.csv` is not the basis for the reported `lower`/`upper` CI columns — those come from the (more accurate) link-scale calculation.
- `car` is loaded but not directly used; `lme4` and `broom.mixed` are loaded a second time, redundantly, before the power-analysis section — both harmless.
