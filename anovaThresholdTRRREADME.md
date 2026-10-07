# Time to Recovery (TTR) Analysis — Threshold Method

R script for estimating how long fish take to return to near-resting oxygen consumption (MO2) after exercise, and comparing that recovery time among four treatment groups.

**Script:** `TTR_threshold_analysis_annotated.R`

---

## Overview

For each fish, the script finds the **time to recovery (TTR)**: the first minute after exercise at which MO2 drops to or below a set percentage of that fish's resting MO2 *and stays there* for the rest of the recovery window. TTR values are then compared among groups with a one-way ANOVA, assumption checks, and Tukey-adjusted post hoc tests, and plotted as a boxplot.

The script analyzes one locomotor mode (walking or swimming) per run.

---

## Requirements

- R (4.0 or later recommended)
- Packages:
  - `tidyverse` — data wrangling and plotting
  - `car` — Levene's test
  - `emmeans` — estimated marginal means and post hoc comparisons

Install with:

```r
install.packages(c("tidyverse", "car", "emmeans"))
```

---

## Input data

`recovery_spreadsheet_1min.csv` — one row per fish per minute of recovery. Required columns:

| Column        | Description                                         |
|---------------|-----------------------------------------------------|
| `Tank_ID`     | Unique fish identifier                              |
| `Mode`        | Locomotor mode: `walking` or `swimming`             |
| `Treatment`   | Rearing environment: `terrestrial` or `aquatic`     |
| `Exercise`    | Exercise treatment: `ramp` or `noRamp`              |
| `minute`      | Minute of recovery (1-min resolution)               |
| `mo2`         | Measured MO2 at that minute                         |
| `resting_mo2` | That fish's resting MO2                             |

Place the CSV in the working directory, or edit `data_file` to point to its location.

---

## Groups

`Treatment` and `Exercise` are combined into four groups:

| Group | Rearing      | Exercise         | Plot color  |
|-------|--------------|------------------|-------------|
| TT    | Terrestrial  | Trained (ramp)   | Dark brown  |
| TU    | Terrestrial  | Untrained (noRamp) | Tan       |
| AT    | Aquatic      | Trained (ramp)   | Dark blue   |
| AU    | Aquatic      | Untrained (noRamp) | Sky blue  |

---

## User settings

Edit these at the top of the script before running:

| Setting         | Default      | Description                                                      |
|-----------------|--------------|------------------------------------------------------------------|
| `data_file`     | `recovery_spreadsheet_1min.csv` | Path to the input CSV                         |
| `selected_mode` | `"walking"`  | Mode to analyze; comment/uncomment to switch to `"swimming"`     |
| `threshold`     | `1.10`       | Recovery cutoff as a multiple of resting MO2 (1.10 = within 10%) |

The recovery window is fixed at **30 minutes** (`filter(minute <= 30)`); change this line to use a different window.

---

## How it works

1. **Load and filter** — reads the CSV, keeps the selected mode, and restricts data to minutes 0–30.
2. **Assign groups** — builds the TT/TU/AT/AU variable.
3. **Set thresholds** — computes each fish's cutoff as `resting_mo2 × threshold`.
4. **Flag recovery** — marks each minute as recovered (`mo2 <= cutoff`) or not.
5. **Calculate TTR** — for each fish, finds the first minute from which every remaining minute is recovered. Fish that never meet this get `TTR = NA`.
6. **Remove non-recovered fish** — drops `NA` rows before analysis.
7. **Transform** — swimming TTR is log-transformed; walking TTR is analyzed raw.
8. **Statistics**
   - One-way ANOVA: `response ~ group`
   - Shapiro–Wilk test on residuals (normality)
   - Levene's test, median-centered (equal variance)
   - Tukey-adjusted pairwise comparisons via `emmeans`
9. **Plot** — boxplot of raw TTR by group with individual fish overlaid.

---

## Outputs

Printed to the console:

- Per-fish TTR table (`NA` = did not recover within the window)
- ANOVA table
- Shapiro–Wilk and Levene's test results
- Estimated marginal means and pairwise contrasts

Plotted:

- Boxplot of TTR by group

Optional: uncomment the `write.csv()` block to save the per-fish TTR table as `TTR_threshold_<mode>.csv`.

---

## Known limitations

- **Non-recovering fish are excluded.** These are the slowest recoverers, so dropping them can bias group means and reduces sample size. A survival analysis treating them as censored at 30 min (e.g., `survival::survdiff`, `coxph`) uses all fish.
- **Edge effect at minute 30.** At the final minute there are no later points to check, so a single below-threshold reading at minute 30 gives `TTR = 30`. Requiring several consecutive recovered minutes would make this more robust.
- **Ceiling effect.** With a 30-min window and 10% threshold, many fish either do not recover or recover only at minutes 28–30, which limits the analysis's ability to detect group differences. Consider a longer window or alternative metrics (e.g., EPOC, MO2 at fixed timepoints).
- **Different scales across modes.** Swimming is log-transformed and walking is not; use the same scale for both if comparing modes.
- **One-way design.** The ANOVA tests the combined 4-level group factor. To test rearing and exercise effects and their interaction separately, use `aov(response ~ Treatment * Exercise, data = TTR_df)`.
- **Harmless plot warning.** `scale_color_manual()` has no mapped color aesthetic. Remove it, or add `aes(color = group)` to `geom_jitter()` to color points.

---

## Usage

1. Set `data_file`, `selected_mode`, and `threshold`.
2. Run the full script in R or RStudio.
3. Switch `selected_mode` and rerun to analyze the other mode.
