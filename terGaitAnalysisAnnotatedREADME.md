# Terrestrial Gait / Kinematics ANOVA Pipeline

This R script (`terGaitAnalysisAnnotated.R`) runs a batch statistical analysis over a set of gait and body-kinematics variables, testing for effects of rearing environment and training on each one. It reads a combined per-trial dataset, cleans and (where needed) transforms each variable, fits a two-way ANOVA, checks assumptions, computes effect sizes and power, runs post-hoc comparisons where warranted, and writes four summary CSVs plus per-variable QQ plots.

## What it does, end to end

For each variable in the analysis list, the script:

1. Selects that variable along with the grouping columns (`group`, `rearedCondition`, `Training`) and drops missing values.
2. Removes within-group outliers using the standard Tukey rule (values outside Q1 − 1.5×IQR to Q3 + 1.5×IQR, computed separately per `group`).
3. Applies the transform specified for that variable — `raw`, `log`, or `logit` — to produce the value actually modeled (`y_used`).
4. Runs Levene's test for equal variances and a two-way ANOVA (`rearedCondition * Training`).
5. Runs a Shapiro-Wilk test on the ANOVA residuals to check normality.
6. Computes eta², Cohen's f, and observed statistical power for each ANOVA term (main effects and interaction).
7. Runs post-hoc pairwise comparisons (via `emmeans`) only where an effect was significant (p < .05) — comparing simple effects within levels of the other factor if the interaction is significant, or comparing each main effect on its own otherwise.
8. Generates a QQ plot of residuals for a visual normality check.

Results from every variable are stacked together and written out at the end, with `rearedCondition` and `Training` relabeled to the friendlier `environment` and `training` in the output tables.

## Input

| File | Description |
|---|---|
| `combinedTerDataSend.csv` | One row per trial/individual. Must contain each variable listed in `vars` (see below) plus the grouping columns `group`, `rearedCondition`, and `Training`. |

Default input path (edit at the top of the script to change):
```
combinedTerDataSend.csv
```

## Variables analyzed

| Variable | Meaning | Transform applied |
|---|---|---|
| `headElevCycleMean_BL` | Mean head elevation per cycle (body lengths) | raw |
| `headElevClipMax_BL` | Max head elevation within a clip (body lengths) | log |
| `curvDorsalCycleMean` | Mean dorsal (body) curvature per cycle | raw |
| `curvDorsalClipMax` | Max dorsal curvature within a clip | raw |
| `mean_stride_length_BL` | Mean stride length (body lengths) | raw |
| `mean_stride_speed_BLs` | Mean stride speed (body lengths/s) | raw |
| `mean_dutyFactor` | Mean duty factor (proportion of stride in stance) | logit |
| `mean_stanceDuration_s` | Mean stance-phase duration (s) | raw |
| `mean_swingDuration_s` | Mean swing-phase duration (s) | log |
| `mean_strideDuration_s` | Mean full stride duration (s) | raw |
| `meanComSpeed_BLs` | Mean center-of-mass speed (body lengths/s) | raw |
| `max_speed_BLs` | Max speed observed (body lengths/s) | raw |
| `median_speed_BLs` | Median speed observed (body lengths/s) | raw |

The `analysis_plan` table at the top of the script is the single place to change a variable's transform. `log` is used for right-skewed, strictly-positive measures; `logit` is used for `mean_dutyFactor` because it's a proportion bounded in [0, 1]. Everything else defaults to `raw`.

## Statistical design

Two-way factorial ANOVA on each (possibly transformed) outcome:

```
y_used ~ rearedCondition * Training
```

- **rearedCondition** — environment the subject was reared in (main effect, labeled `environment` in output)
- **Training** — training condition/group (main effect, labeled `training` in output)
- **rearedCondition:Training** — interaction between the two

Assumption checks reported per variable: Levene's test (variance homogeneity) and Shapiro-Wilk on residuals (normality).

Post-hoc comparisons only run where an effect is significant (α = .05):
- Interaction significant → simple effects of `rearedCondition` within each `Training` level, and of `Training` within each `rearedCondition` level.
- Interaction not significant → each significant main effect compared on its own via `emmeans`/`pairs()`.

Effect size and power, per ANOVA term:
- **eta²** = SS(effect) / SS(total)
- **Cohen's f** = √(eta² / (1 − eta²))
- **observed power** via `pwr::pwr.f2.test()` at α = .05

## Output

Written to the same directory as the input file:

| File | Contents |
|---|---|
| `terStats.csv` | ANOVA table (effect, df, sum of squares, F, p) for every variable |
| `terPosthoc.csv` | Pairwise post-hoc comparisons, only for effects that were significant |
| `terDiagnostics.csv` | Shapiro-Wilk and Levene results (and the transform used) per variable |
| `terPower.csv` | eta², Cohen's f, and observed power for each ANOVA term per variable |

Default output directory (edit at the bottom of the script to change):
```
dataOutputPolyp
```

A QQ plot of residuals is also generated for each variable in the analysis loop (printed to the active graphics device — wrap the script's loop in `pdf(...)` / `dev.off()` if you want these saved to a file rather than shown interactively).

## Dependencies

R packages: `readr`, `dplyr`, `car`, `rstatix`, `broom`, `emmeans`, `nlme`, `ggplot2`, `pwr`.

Install any missing packages with:
```r
install.packages(c("readr", "dplyr", "car", "rstatix", "broom",
                    "emmeans", "nlme", "ggplot2", "pwr"))
```

## Running it

1. Update the input path near the top of the script (and the `out_dir` path near the bottom) if your file locations differ.
2. Run the script top to bottom in R or RStudio (`source("terGaitAnova_annotated.R")`).
3. Check console output for the per-variable ANOVA summary printed by `print(anova_results)`, and inspect the four CSVs written to `out_dir`.

## Notes

- The `logit` transform clips values to `(1e-6, 1 − 1e-6)` before transforming, so exact 0s or 1s in `mean_dutyFactor` won't produce infinite values.
