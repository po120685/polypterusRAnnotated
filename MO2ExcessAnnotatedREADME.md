# Total Energy Analysis — Swimming vs. Walking

This README documents `MO2ExcessAnnotated.R`, which analyzes total
energy expenditure (excess post-exercise oxygen consumption, converted to
energy units) for fish tested under swimming and walking conditions, across
two experimental factors: **Treatment** (aquatic vs. terrestrial rearing) and
**Exercise** (ramped exercise vs. no ramp).

## What the script does

For each locomotion **Mode** (swimming, walking) analyzed separately, the
script:

1. Runs a classic 2-way ANOVA (`Treatment * Exercise`) on `total_energy_J`.
2. Runs a permutation-based 2-way ANOVA (`permuco::aovperm`, 5000
   permutations) as a robustness check that doesn't rely on normality.
3. Checks ANOVA assumptions: Shapiro-Wilk test for normality of residuals,
   Levene's test for equal variance, and the standard `plot.lm` diagnostic
   panel.
4. Runs post-hoc pairwise comparisons across all four Treatment × Exercise
   groups via `emmeans`.
5. Makes boxplots (with jittered raw points) of total energy by group, one
   per Mode, using a fixed color scheme and group order.
6. Saves the permutation ANOVA result tables to CSV.
7. Runs a simulation-based power analysis for the `Treatment:Exercise`
   interaction term specifically — sweeping across a range of hypothetical
   effect sizes (0.25x–3x the observed effect) to estimate how much power
   the current sample size has to detect that interaction, and at what
   effect size 80% power would be reached.

## Requirements

R packages: `tidyverse`, `car`, `emmeans`, `permuco`, `svglite`, `future`,
`future.apply`. Install with:

```r
install.packages(c("tidyverse", "car", "emmeans", "permuco", "svglite",
                    "future", "future.apply"))
```

## Input file

`MO2TTRSumPerFish.csv` — one row per fish per Mode trial. Columns:

| Column | Meaning |
|---|---|
| `fishID` | Unique fish identifier, `Tank_ID_Fish_no` (e.g. `A-04_26`). |
| `Tank_ID` | Tank/housing unit the fish was reared in. |
| `Fish_no` | Sequential number assigned to the individual fish. |
| `Mode` | Locomotion trial type: `swimming` or `walking`. |
| `Treatment` | Rearing/acclimation condition: `aquatic` or `terrestrial`. |
| `Exercise` | Whether the fish underwent the ramped exercise protocol (`ramp`) or a resting/control protocol (`noRamp`). |
| `group` | 2-letter combined code: 1st letter = Treatment (A/T), 2nd letter = Exercise (T = trained/ramp, U = untrained/noRamp). So `AT`, `AU`, `TT`, `TU`. This is what the boxplots color by. |
| `TTR` | "Time To Recovery": recovery-period duration (hours) used as the integration window for the excess totals below. Fixed per Mode × Treatment × Exercise cell rather than measured per fish — appears to be a protocol parameter, not a per-fish outcome. |
| `total_MO2_excess` | Total oxygen consumption in excess of resting baseline, integrated over the TTR window (i.e. EPOC — excess post-exercise oxygen consumption). |
| `total_energy_cal` | `total_MO2_excess` converted to calories using an oxycaloric equivalent of 3.25 cal per unit MO2 excess (consistent with RQ ≈ 0.8; confirmed exactly against the data). |
| `total_energy_J` | `total_energy_cal` converted to joules (× 4.184 J/cal; confirmed exactly against the data). **This is the outcome variable modeled throughout the script.** |


## Outputs

All written to `out_dir` (`dataOutputPolyp`):

- `swimmingEnergyExcessBoxplot.svg`, `walkingEnergyExcessBoxplot.svg` —
  boxplots of `total_energy_J` by group, one per Mode.
- `swimmingEnergyPermANOVA.csv`, `walkingEnergyPermANOVA.csv` —
  permutation ANOVA result tables (F-statistics, resampled p-values) per
  Mode.
- `powerCurveInteraction.csv` — power, 95% CI, and number of valid
  simulation replicates for the `Treatment:Exercise` interaction, at each
  tested effect-size multiplier, for both modes.
- `powerCurve_Interaction_Po.svg` — power curve plot (power vs. effect size,
  one line per Mode, with a dashed reference line at 80% power).

Console-only output (not saved to file): `summary()` of each ANOVA model,
Shapiro-Wilk and Levene test results, the `plot.lm` diagnostic panels, and
the `emmeans` pairwise comparison tables. If you want these preserved,
consider redirecting console output (`sink()`) or exporting the emmeans
objects to CSV as well.

## How the power analysis works

`run_power_curve()` fits a linear model to the *real* data to extract
realistic Treatment/Exercise/interaction coefficients and residual noise,
then for each effect-size multiplier (0.25x–3x) it:

1. Simulates 200 synthetic datasets by scaling those coefficients and adding
   fresh random noise (same residual SD as observed).
2. Re-runs the permutation ANOVA (299 permutations, for speed) on each
   simulated dataset.
3. Records the fraction of simulations where the `Treatment:Exercise`
   p-value was significant (< 0.05) — that fraction is the power estimate,
   with a binomial 95% CI.

Simulation runs in parallel via `future`/`future.apply` (`multisession`,
safe on RStudio/macOS). This step is the most computationally expensive part
of the script (200 simulations × 7 effect sizes × 2 modes = 2800 permutation
ANOVA fits) and may take a while depending on core count.