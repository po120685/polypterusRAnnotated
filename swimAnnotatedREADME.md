# Swim Master Analysis Pipeline

This R script (`swimAnnotated.R`) runs a mixed-model / permutation-ANOVA / power-analysis pipeline across 15 swimming and kinematic variables (wave, tail/fin, effort, and hydrodynamic metrics), reading from and writing to the `swimSend` data folder.

## What it does, end to end

For each of 15 outcome variables, using a manually assigned method (`LMM`, `logLMM`, or `perm`):

1. **LMM** — fits `y ~ (environment + training + flowSpeed)^2 + (1 | fishID)` with `lme4`/`lmerTest`, prints the model summary and ANOVA table, runs post-hoc comparisons on any significant two-way interaction (all three interactions handled symmetrically), records residual-normality and heteroscedasticity diagnostics, and stores the fixed-effect coefficients.
2. **logLMM** — identical to LMM, but fit on `log(y)` (for positive, right-skewed variables); rows are filtered to `y > 0` first.
3. **perm** — fits a permutation ANOVA (`permuco::aovperm`) on `y ~ (environment + training + flowSpeed)^2`, with permutations clustered by `fishID`; additionally runs a post-hoc test if `environment:training` is significant (see below).

For every variable, the script also runs a **power analysis** for each fixed-effect term:
- LMM / logLMM → simulation-based power via `simr::powerSim()`, using the model's own fitted effect size.
- perm → a custom simulate-and-count procedure: inject a known effect (1×, 2×, and 3× the LMM-estimated coefficient) into synthetic data, rerun the permutation ANOVA 1,000 times per effect size, and report the fraction of runs with p < .05, with a binomial exact confidence interval.

Results from all 15 variables are stacked into four tables and written out as CSVs.

## Input

| File | Description |
|---|---|
| `combinedSwimDataSend.csv` | One row per trial/timepoint, with repeated rows per `fishID`. Must contain each variable in `analysis_plan` plus `fishID`, `environment` (0/1), `training` (0/1), and `flowSpeed`. |

Default input path (edit near the top of the script to change):
```
combinedSwimDataSend.csv
```

`environment` and `training` are recoded from 0/1 into readable factor labels: `0 = Aquatic`, `1 = Terrestrial` for environment; `0 = Untrained`, `1 = Trained` for training.

## Variables analyzed

| Variable | Meaning | Method |
|---|---|---|
| `strideLength_BL` | Stride length (body lengths) | perm |
| `meanWaveSpeed_BLs` | Mean body-wave propagation speed (body lengths/s) | LMM |
| `meanWavelength_BL` | Mean body-wave wavelength (body lengths) | perm |
| `meanWavePeriod_s` | Mean body-wave period (s) | logLMM |
| `meanWaveFreq_Hz` | Mean body-wave frequency (Hz) | LMM |
| `meanTailFreq_Hz` | Mean tail-beat frequency (Hz) | LMM |
| `meanFinFreq_Hz` | Mean fin-beat frequency (Hz) | LMM |
| `meanTailAmp_BL` | Mean tail-beat amplitude (body lengths) | logLMM |
| `meanFinAmp_BL` | Mean fin-beat amplitude (body lengths) | perm |
| `finEffort` | Fin-based locomotor effort | perm |
| `tailEffort` | Tail-based locomotor effort | logLMM |
| `totalEffort` | Combined locomotor effort | perm |
| `Re` | Reynolds number (inertial vs. viscous forces) | perm |
| `strouhal` | Strouhal number (dimensionless oscillation-frequency parameter, common in swimming biomechanics) | logLMM |
| `angle_deg` | Body angle (degrees) | perm |

## Statistical design

All three methods test the same fixed-effects structure:

```
y ~ (environment + training + flowSpeed)^2 [+ (1 | fishID) for LMM/logLMM]
```

The `^2` shorthand expands to all three main effects plus all three pairwise interactions (`environment:training`, `environment:flowSpeed`, `training:flowSpeed`). `(1 | fishID)` is a per-fish random intercept, accounting for the fact that each fish contributes multiple rows.

### Post-hoc comparisons

**LMM/logLMM** — runs per two-way interaction, only when that interaction is significant (p < .05), and treats all three interactions the same way:
- `environment:training` significant → pairwise comparison of estimated marginal means across the four environment×training cells (Tukey-adjusted).
- `environment:flowSpeed` or `training:flowSpeed` significant → comparison of the flowSpeed *slope* between environment (or training) levels, via `emmeans::emtrends()`, since flowSpeed is continuous.

**perm** — only follows up `environment:training`; `environment:flowSpeed` and `training:flowSpeed` are never given a post-hoc step in this branch, even if significant (the p-values for those two are computed but not checked). When `environment:training` *is* significant, the follow-up refits the data as an ordinary linear model (`lm(y ~ environment * training + flowSpeed)`, no random effect) and runs Tukey-adjusted pairwise comparisons of environment within each training level via `emmeans`. In other words: the omnibus test is permutation-based, but its post-hoc follow-up is parametric — a practical compromise, since the permutation ANOVA itself has no direct pairwise-comparison analog, but worth keeping in mind when interpreting those specific p-values.

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
| `swimStats.csv` | Model results for every variable — LMM/logLMM fixed-effect coefficients (estimate, SE, df, statistic, p-value) or permutation-ANOVA statistics (F, p-value; estimate/SE/df left NA) |
| `swimDiagnostics.csv` | Shapiro-Wilk and Breusch-Pagan results per variable (NA for `perm`-method variables) |
| `swimPosthoc.csv` | Pairwise/slope post-hoc comparisons — all significant interactions for LMM/logLMM; `environment:training` only for `perm` |
| `swimPower.csv` | Power, 95% CI, and effect estimate per variable × fixed-effect term (and, for `perm` variables, per effect-size multiplier) |

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
2. Run the script top to bottom in R or RStudio (`source("swimMasterAnalysis_annotated.R")`). This is computationally heavy — with 15 variables, several of which use `perm` (5,000 permutations for the main test, plus 1,000 simulated datasets × up to 6 terms × 3 effect sizes for power, each running its own 5,000-permutation test), expect a long runtime.
3. Check console output for per-variable model summaries and post-hoc results as they print, and inspect the four CSVs written to the output directory.

