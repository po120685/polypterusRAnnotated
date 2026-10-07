# Time to Recovery (TTR) Analysis — Two-Way ANOVA

R script for estimating how long fish take to return to near-resting oxygen consumption (MO2) after exercise, and testing whether recovery time depends on rearing environment, exercise training, or their interaction.

**Script:** `TTR_threshold_twoway_annotated.R`

---

## Overview

For each fish, the script finds the **time to recovery (TTR)**: the first minute after exercise at which MO2 drops to or below a set percentage of that fish's resting MO2 *and stays there* for the rest of the recovery window.

TTR is then analyzed with a **2 × 2 two-way ANOVA**:

- **Treatment** (rearing environment): terrestrial vs aquatic
- **Exercise** (training): ramp vs noRamp
- **Treatment × Exercise** interaction

The script analyzes one locomotor mode (walking or swimming) per run.

---

## Requirements

- R (4.0 or later recommended)
- Packages:
  - `tidyverse` — data wrangling and plotting
  - `car` — Type III ANOVA (`Anova()`) and Levene's test
  - `emmeans` — estimated marginal means and post hoc comparisons

```r
install.packages(c("tidyverse", "car", "emmeans"))
```

---

## Input data

`recovery_spreadsheet_1min.csv` — one row per fish per minute of recovery.

| Column        | Description                                       |
|---------------|---------------------------------------------------|
| `Tank_ID`     | Unique fish identifier                            |
| `Mode`        | Locomotor mode: `walking` or `swimming`           |
| `Treatment`   | Rearing environment: `terrestrial` or `aquatic`   |
| `Exercise`    | Exercise treatment: `ramp` or `noRamp`            |
| `minute`      | Minute of recovery (1-min resolution)             |
| `mo2`         | Measured MO2 at that minute                       |
| `resting_mo2` | That fish's resting MO2                           |

Place the CSV in the working directory, or edit `data_file` to point to it.

---

## Groups

| Group | Rearing     | Exercise           | Plot color |
|-------|-------------|--------------------|------------|
| TT    | Terrestrial | Trained (ramp)     | Dark brown |
| TU    | Terrestrial | Untrained (noRamp) | Tan        |
| AT    | Aquatic     | Trained (ramp)     | Dark blue  |
| AU    | Aquatic     | Untrained (noRamp) | Sky blue   |

The ANOVA uses `Treatment` and `Exercise` directly; the `group` label is used for plots and tables.

---

## User settings

| Setting         | Default                          | Description                                              |
|-----------------|----------------------------------|----------------------------------------------------------|
| `data_file`     | `recovery_spreadsheet_1min.csv`  | Path to the input CSV                                    |
| `selected_mode` | `"walking"`                      | Mode to analyze; comment/uncomment to switch to `"swimming"` |
| `threshold`     | `1.10`                           | Recovery cutoff as a multiple of resting MO2 (1.10 = within 10%) |
| `max_minute`    | `30`                             | Length of the recovery window (minutes)                  |
| `use_log`       | `TRUE` for swimming, `FALSE` for walking | Whether to analyze log(TTR) or raw TTR           |

**Recommended:** set `use_log <- FALSE` so both modes are analyzed on the same (raw) scale. TTR values span only ~28–30 min, so a log transform has almost no effect, and residuals are approximately normal without it. If you keep the log for swimming, state the reason in your methods.

---

## How it works

1. **Load and filter** — reads the CSV, keeps the selected mode, restricts to minutes 0–`max_minute`.
2. **Set factors and groups** — fixes level order for `Treatment` and `Exercise`, builds TT/TU/AT/AU labels.
3. **Flag recovery** — each minute is marked recovered if `mo2 <= resting_mo2 × threshold`.
4. **Calculate TTR** — first minute from which all remaining minutes are recovered; `NA` if never.
5. **Summarize recovery** — counts recovered and non-recovered fish per group.
6. **Remove non-recovered fish** and print cell sizes entering the ANOVA.
7. **Choose the response scale** — log or raw, per `use_log`.
8. **Two-way ANOVA** — `aov(TTR ~ Treatment * Exercise)`, or `aov(log(TTR) ~ Treatment * Exercise)` when `use_log = TRUE`, tested with `car::Anova(type = 3)` using sum-to-zero contrasts. The log is written inside the model formula so `emmeans` can back-transform estimates to minutes.
9. **Assumption checks** — Shapiro–Wilk and Q-Q plot of residuals; Levene's test across the four cells.
10. **Post hoc tests** (`emmeans`, Tukey-adjusted):
    - Main effects of Treatment and Exercise
    - Simple effects (Exercise within each Treatment, and vice versa)
    - All pairwise comparisons among the four cells
11. **Plots** — boxplot of TTR by group and an interaction plot of cell means ± SE.

---

## Outputs

Printed to the console:

- Per-fish TTR table (`NA` = did not recover within the window)
- Recovery summary (recovered vs not, per group)
- Cell sizes used in the ANOVA
- Two-way ANOVA table (Type III)
- Shapiro–Wilk and Levene's test results
- Estimated marginal means and pairwise contrasts

Plotted:

- Q-Q plot of residuals
- Boxplot of TTR by group
- Treatment × Exercise interaction plot

Optional: uncomment the `write.csv()` block to save the per-fish TTR table (including non-recovered fish) as `TTR_threshold_<mode>.csv`.

---

## Interpreting the results

- **ANOVA table** — read the `Treatment`, `Exercise`, and `Treatment:Exercise` rows; ignore `(Intercept)`. Report each as F(1, df_residual) and p.
- **If the interaction is significant** — interpret the simple effects, not the main effects.
- **If the interaction is not significant** — interpret the main effects.
- **Log-transformed analyses** — emmeans are back-transformed to minutes (geometric means), and contrasts are reported as ratios.
- **Identical F-values for Exercise and the interaction** — this can happen when the exercise effect is exactly zero in one rearing group (e.g., both terrestrial groups have the same TTR). It is a property of the data, not an error.
- **Reporting format** — e.g., "TTR was not affected by rearing environment (two-way ANOVA, F(1,6) = 0.34, p = 0.58), exercise treatment (F(1,6) = 0.93, p = 0.37), or their interaction (F(1,6) = 0.04, p = 0.85)."
- **Interaction plot** — parallel lines suggest no interaction; crossing or diverging lines suggest the exercise effect differs by rearing environment.

---

## Why Type III sums of squares

Removing non-recovered fish leaves unequal group sizes. With unequal groups, base R's `summary(aov())` uses sequential (Type I) sums of squares, so results change with term order. Type III tests each effect adjusted for all others and requires sum-to-zero contrasts, which the script sets with `options(contrasts = c("contr.sum", "contr.poly"))`.

If the interaction is not significant, Type II (`Anova(model, type = 2)`) is a reasonable alternative for testing main effects.

---

## Usage

1. Set `data_file`, `selected_mode`, `threshold`, `max_minute`, and `use_log`.
2. Run the full script in R or RStudio.
3. Switch `selected_mode` and rerun for the other mode.

---

## Suggested methods text

> Time to recovery (TTR) was defined as the first minute at which MO2 fell to and remained at or below 110% of resting MO2 within a 30-minute recovery period. Fish that did not recover were excluded. For each locomotor mode, TTR was analyzed with a two-way ANOVA with rearing environment, exercise treatment, and their interaction as fixed effects, using Type III sums of squares to account for unequal sample sizes. Pairwise differences were assessed with Tukey-adjusted comparisons of estimated marginal means. Normality of residuals was evaluated with a Shapiro–Wilk test and homogeneity of variance with Levene's test. Analyses were performed in R using the car and emmeans packages.
