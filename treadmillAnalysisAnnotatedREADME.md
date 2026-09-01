# Treadmill Master Analysis Pipeline

This R script (`treadmillMasterAnalysisAnnotated.R`) is the mixed-model counterpart to the terrestrial-gait ANOVA pipeline: it analyzes 8 treadmill/kinematic outcome variables that have **repeated measurements per fish**, so each variable is modeled with `fishID` as a random effect (or, for two variables, with a permutation test that clusters permutations by fish) rather than a simple one-row-per-subject ANOVA.

## What it does, end to end

For each of 8 outcome variables, using a manually assigned method (`LMM`, `logLMM`, or `perm`):

1. **LMM** — fits `y ~ (environment + training + flowSpeed)^2 + (1 | fishID)` with `lme4`/`lmerTest`, prints the model summary and ANOVA table, runs post-hoc comparisons on any significant two-way interaction, records residual-normality and heteroscedasticity diagnostics, and stores the fixed-effect coefficients.
2. **logLMM** — identical to LMM, but fit on `log(y)` (for positive, right-skewed variables); rows are filtered to `y > 0` first.
3. **perm** — fits a permutation ANOVA (`permuco::aovperm`) on `y ~ (environment + training + flowSpeed)^2`, with permutations clustered by `fishID` (used when a variable doesn't fit the LMM's normality/variance assumptions well).

For every variable, the script also runs a **power analysis** for each fixed-effect term:
- LMM / logLMM → simulation-based power via `simr::powerSim()`, using the model's own fitted effect size.
- perm → a custom simulate-and-count procedure: inject a known effect (1×, 2×, and 3× the LMM-estimated coefficient) into synthetic data, rerun the permutation ANOVA 1,000 times per effect size, and report the fraction of runs with p < .05, with a binomial exact confidence interval.

Results from all 8 variables are stacked into four tables and written out as CSVs.

## Input

| File | Description |
|---|---|
| `combinedWalkDataSend.csv` | One row per trial/timepoint, with repeated rows per `fishID`. Must contain each variable in `analysis_plan` plus `fishID`, `environment` (0/1), `training` (0/1), and `flowSpeed`. |

Default input path (edit near the top of the script to change):
```
combinedWalkDataSend.csv
```

`environment` and `training` are recoded from 0/1 into readable factor labels: `0 = Aquatic`, `1 = Terrestrial` for environment; `0 = Untrained`, `1 = Trained` for training.

## Variables analyzed

| Variable | Meaning | Method |
|---|---|---|
| `strideLength_BL` | Stride length (body lengths) | perm |
| `finEffort` | Fin-based locomotor effort | LMM |
| `tailEffort` | Tail-based locomotor effort | perm |
| `totalEffort` | Combined locomotor effort | perm |
| `tailMeanAmp_BL` | Mean tail oscillation amplitude (body lengths) | LMM |
| `finMeanAmp_BL` | Mean fin oscillation amplitude (body lengths) | perm |
| `tailMeanFreq_Hz` | Mean tail oscillation frequency (Hz) | logLMM |
| `finMeanFreq_Hz` | Mean fin oscillation frequency (Hz) | logLMM |

The `analysis_plan` table near the top of the script is the single place to change a variable's method — these choices were made manually (presumably from prior diagnostic checks), not derived automatically inside the script itself.

## Statistical design

All three methods test the same fixed-effects structure:

```
y ~ (environment + training + flowSpeed)^2 [+ (1 | fishID) for LMM/logLMM]
```

The `^2` shorthand expands to all three main effects plus all three pairwise interactions (`environment:training`, `environment:flowSpeed`, `training:flowSpeed`), but **not** the three-way interaction. `(1 | fishID)` is a per-fish random intercept, accounting for the fact that each fish contributes multiple rows.

**Post-hoc comparisons** (LMM/logLMM only) run per two-way interaction, only when that interaction is significant (p < .05):
- `environment:training` significant → pairwise comparison of estimated marginal means across the four environment×training cells (Tukey-adjusted).
- `environment:flowSpeed` or `training:flowSpeed` significant → comparison of the flowSpeed *slope* between environment (or training) levels, via `emmeans::emtrends()`, since flowSpeed is continuous and a simple mean comparison doesn't apply.

**Diagnostics** (LMM/logLMM only): Shapiro-Wilk normality and a Breusch-Pagan-style heteroscedasticity test on model residuals, both via the `performance` package. `perm`-method variables get `NA` diagnostics rows, since a permutation test makes no distributional assumptions to check.

**Effect-name standardization**: `lmer` names dummy-coded terms after their factor level (e.g. `environmentTerrestrial`), while `aovperm` reports plain factor names (e.g. `environment`). The script rewrites LMM/logLMM term names to match the permutation ANOVA's naming (`environment`, `training`, `environment:training`, …) so results from all three methods can be filtered/compared by the same `effect` label in the output tables.

## Power analysis

- **LMM / logLMM**: `simr::powerSim()` resimulates the response from the fitted model 1,000 times per fixed-effect term and refits, reporting the proportion of simulations reaching significance (parsed from `powerSim`'s printed summary, since it isn't exposed as a clean numeric field).
- **perm**: since `permuco` has no built-in power function, the script does its own simulate-and-count loop per term:
  1. Fits a helper LMM (not saved to the main results) just to get a realistic effect-size estimate for each term.
  2. Tests injected effect sizes of 1×, 2×, and 3× that estimate.
  3. At each size, builds 1,000 synthetic datasets (`simulate_data()` — real baseline values + injected effect + per-fish random intercept + residual noise scaled to the real data), reruns the permutation ANOVA on each (`run_perm_test()`), and computes power as the fraction with p < .05, with a binomial exact CI (`binom.test()`).

A fixed `set.seed(123)` is used throughout so permutation tests and power simulations are reproducible across runs.

## Output

Written to the same directory as the input file:

| File | Contents |
|---|---|
| `treadmillStats.csv` | Model results for every variable — LMM/logLMM fixed-effect coefficients (estimate, SE, df, statistic, p-value) or permutation-ANOVA statistics (F, p-value; estimate/SE/df left NA) |
| `treadmillDiagnostics.csv` | Shapiro-Wilk and Breusch-Pagan results per variable (NA for `perm`-method variables) |
| `treadmillPosthoc.csv` | Pairwise/slope post-hoc comparisons, only for interactions that were significant (LMM/logLMM only) |
| `treadmillPower.csv` | Power, 95% CI, and effect estimate per variable × fixed-effect term (and, for `perm` variables, per effect-size multiplier) |

Default output directory (edit near the bottom of the script to change):
```
dataOutputPolyp
```

## Dependencies

R packages: `readr`, `dplyr`, `lme4`, `lmerTest`, `emmeans`, `performance`, `permuco`, `simr`, `tidyr`.

Install any missing packages with:
```r
install.packages(c("readr", "dplyr", "lme4", "lmerTest", "emmeans",
                    "performance", "permuco", "simr", "tidyr"))
```

## Running it

1. Update the input path near the top of the script, and the output paths near the bottom, if your file locations differ.
2. Run the script top to bottom in R or RStudio (`source("treadmillMasterAnalysis_annotated.R")`). This is computationally heavy — the `perm` variables alone run 5,000 permutations for the main test plus 1,000 simulated datasets × up to 6 terms × 3 effect sizes for power (each running its own 5,000-permutation test), so expect a long runtime.
3. Check console output for per-variable model summaries and post-hoc results as they print, and inspect the four CSVs written to the output directory.

## Notes / things to double-check when reading results

- `perm`-method power estimates depend on effect sizes borrowed from a helper LMM fit on the raw (untransformed, non-permutation) data — they are not derived from the permutation model itself, since `aovperm` doesn't produce comparable coefficient estimates.
- `logLMM` power analysis is refit on a *new* log-transformed column (`target_var`) inside the main loop, separately from the log-transformed fit inside `run_log_lmm()` — both should agree, but they are two separate model fits.
- In `run_lmm()`, three of the "standardize effect names" `gsub()` calls replace a pattern with itself (e.g. `gsub("environment:training", "environment:training", ...)`) — these are no-ops left in as an explicit marker that those names are already in the target format; they don't affect the output.
- Post-hoc and diagnostics tables only contain rows for variables/effects where the relevant test applies — a `perm`-method variable won't appear in `treadmillPosthoc.csv` or with real values in `treadmillDiagnostics.csv`.
