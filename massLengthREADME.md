# Polypterus body size analyses

R code testing whether body mass and length differ between rearing environments
(aquatic vs. terrestrial) and training treatments (ramp vs. no ramp) in
*Polypterus* used in the endurance and kinematics experiments.

## Files

| File | Description |
|---|---|
| `polypterus_size_anovas.R` | Analysis script |
| `fishMassLength.csv` | Input data: one row per fish per trial (84 trials, 24 fish) |

The script can also read the original workbook
(`polypterus_all_experiments_2.xlsx`, sheet "All experiments"). Change the
`file <-` line near the top of the script to switch.

## Data

| Column | Units | Definition |
|---|---|---|
| Experiment type | — | Endurance or Kinematics |
| Test | — | Swim, Walk or Full terrestrial (kinematics only, no speed control) |
| Fish ID | — | Individual identifier; aquatic-reared fish have alphanumeric codes, terrestrial-reared fish colour-tag combinations |
| Environment | — | Rearing environment: Aquatic or Terrestrial |
| Ramp | — | Training: R = trained with a ramp, NR = no ramp (untrained) |
| Mass (g) | g | Body mass at the time of the trial |
| Length (cm) | cm | Body length at the time of the trial |
| Body condition K | — | Fulton's condition factor, K = 100 × M / L³ (not used in this script) |

## Requirements

- R (tested with R 4.3.3)
- Packages: `car` (tested with 3.1.2) and `readxl` (tested with 1.4.3; only
  needed if reading the .xlsx)

```r
install.packages(c("car", "readxl"))
```

## How to run

1. Put `polypterus_size_anovas.R` and `fishMassLength.csv` in the same folder.
2. Set that folder as the R working directory, e.g. `setwd("path/to/folder")`.
3. Run the whole script: `source("polypterus_size_anovas.R")`.

## What the script does

The same fish were used in more than one test, so each experiment × test
combination is analysed separately. Each fish therefore appears only once in
each analysis. There are five subsets:

| Experiment | Test | n fish |
|---|---|---|
| Endurance | Swim | 16 |
| Endurance | Walk | 20 |
| Kinematics | Swim | 12 |
| Kinematics | Walk | 12 |
| Kinematics | Full terrestrial | 24 |

For each subset and each response (mass, length), the script fits a two-way
ANOVA:

```
response ~ Environment * Ramp
```

using Type II sums of squares (`car::Anova`), because group sizes are unequal.

It then applies a decision rule for model assumptions (α = 0.05):

1. Fit the model on raw data. Test normality of residuals (Shapiro–Wilk) and
   homogeneity of variance across the four Environment × Ramp groups (Levene's
   test, median-centred / Brown–Forsythe, the `car::leveneTest` default).
2. If either test fails, refit on natural-log-transformed data.
3. If the log-transformed data still fail either test, use an ANOVA with
   HC3 heteroscedasticity-consistent standard errors
   (`Anova(..., white.adjust = "hc3")`).

The model reached at the last step is the one reported.

## Output

Printed to the console, for each subset:

- the number of fish per Environment × Ramp group
- Shapiro–Wilk and Levene p-values at each step of the decision rule
- the ANOVA table for the reported model

Written to the working directory:

| File | Contents |
|---|---|
| `size_anova_results.csv` | One row per term (Environment, Ramp, Environment:Ramp) per subset and response: model used, df, F, p, assumption-test p-values and a `*` where p < 0.05 |
| `size_group_means.csv` | Mean, SD and n of mass and length for each experiment × test × Environment × Ramp group |

## Results with the current data

| Subset | Mass | Length |
|---|---|---|
| Endurance swim | ln-transformed; no significant effects (Environment p = 0.052) | No significant effects |
| Endurance walk | ln-transformed + HC3; no significant effects | No significant effects |
| Kinematics swim | Environment p = 0.026 | Environment p = 0.026 |
| Kinematics walk | Environment p < 0.001; Ramp p = 0.002 | Environment p = 0.012; Ramp p = 0.015 |
| Kinematics full terrestrial | Environment p = 0.039 | No significant effects |

No Environment × Ramp interaction was significant.

## Notes

- Kinematics swim and walk have only 3 fish per group, so these tests have low
  power; a non-significant result does not show that groups were size-matched.
- To use the classic mean-centred Levene's test instead, add `center = mean` to
  the `leveneTest()` call in `fit_anova()`.
- The script stops with an error if a fish appears more than once within a
  subset (`stopifnot(!any(duplicated(sub$Fish)))`).

## Reference

Fox, J. & Weisberg, S. (2019). *An R Companion to Applied Regression*, 3rd ed.
Sage, Thousand Oaks, CA.
