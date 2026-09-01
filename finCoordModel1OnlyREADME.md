# Fin-Use Analysis Pipeline

This R script (`finCoordModel1Only.R`) fits a Bayesian multinomial mixed model of fin-use behavior in the treadmill dataset: `finUse` (a 3-level outcome — None / Asymmetrical / Symmetrical) modeled as a function of environment, training, and flow speed, with a per-fish random intercept.

## What it does, end to end

1. Fits `finUse ~ (environment + training + flowSpeed)^2 + (1 | fishID)` as a categorical (multinomial logistic) model via `brms`/Stan (4 chains × 4,000 iterations), with weakly informative priors and tuned sampler settings to handle a model this complex.
2. Saves the fitted model object to disk (`model.rds`) and runs standard Bayesian diagnostics: `summary()`, trace/density plots, and a posterior predictive check (`pp_check()`).
3. Computes population-level predicted probabilities of each fin-use category for every row, saves them (`finCoordinationModel1.csv`), and aggregates them to one row per fish × environment × training × flowSpeed cell (`finCoordinationModel2.csv`).
4. Plots predicted probability vs. flow speed, faceted by environment × training, one line per category.
5. Runs a large simulation-based power analysis for every fixed effect in the model (see below).

**Naming note:** `finCoordinationModel1.csv` and `finCoordinationModel2.csv` are **both outputs of this same model** — the "Model1"/"Model2" in those filenames refers to two output *tables* (row-level predictions vs. fish-aggregated predictions), not to two different statistical models. There's only one model in this script.

## Input

| File | Description |
|---|---|
| `combinedWalkDataSend.csv` | One row per trial/timepoint, repeated per `fishID`. Must contain `finUse` (a "None"/"Asymmetrical"/"Symmetrical" column) and `base_id`. |

Default input path (edit near the top of the script to change):
```
/Users/theopo/Documents/MATLAB/polyp/dataOutputPolyp/valeSend/walkSend/combinedWalkDataSend.csv
```
`environment` and `training` are recoded from 0/1 into readable factor labels (`Aquatic`/`Terrestrial`, `Untrained`/`Trained`). `finUse` is set up with `"None"` as the first (reference) factor level, which is why coefficient names throughout use `muAsymmetrical_...` / `muSymmetrical_...` prefixes (brms's categorical family reports each non-reference category relative to the reference).

## Statistical design

```
finUse ~ (environment + training + flowSpeed)^2 + (1 | fishID)
```

expanding to all three main effects plus the three pairwise interactions (not the three-way interaction), with a per-fish random intercept.

- `prior = normal(0, 2)` on all fixed-effect coefficients for both non-reference categories — a weakly informative prior to prevent unbounded estimates from near-perfect separation (a common issue in categorical/logistic models with sparse categories).
- `control = list(adapt_delta = 0.99, max_treedepth = 15)` — standard Stan/brms sampler tuning to reduce divergent transitions and treedepth warnings in a model this complex.
- `chains = 4, cores = 4, iter = 4000` for the main reported fit.
- Diagnostics: `summary(model)` (posterior estimates + Rhat/ESS), `plot(model)` (trace/density plots), `pp_check(model)` (posterior predictive check).

## Power analysis

Worth understanding what kind of power analysis this is: it's a *retrospective/self-consistency* question — "if the world really looks exactly like what the fitted model's posterior says, and I collected a brand-new dataset of the same size and refit, how often would each effect's 95% credible interval exclude zero?" It only tests detectability at the model's own observed (fitted) effect size, not a swept range of hypothetical effect sizes.

Mechanically, for each of `nsim` (1,000) iterations:
1. Draw one full posterior-predictive replicate of `finUse` for every row via `posterior_predict(model)`, and pick one posterior draw's worth of simulated categories at random.
2. Swap the real `finUse` column for that simulated one (predictors and `fishID` unchanged).
3. Refit the multinomial model on the simulated dataset, reusing the already-compiled Stan model (`recompile = FALSE`) but with lighter sampler settings than the main fit (2 chains × 2,000 iterations vs. 4 × 4,000), with a fresh seed each time.
4. For every main effect, interaction, and one "terrestrial simple slope" (the effective flowSpeed slope specifically within Terrestrial-reared fish — main effect + interaction coefficient, summed posterior-draw-by-draw), check whether the 95% credible interval excludes zero in this refit.
5. Write a checkpoint CSV after every *successful* iteration (a failed refit skips the checkpoint write for that iteration via `next`, so the file reflects however many simulations most recently succeeded).

The final `power_summary` table reports, per effect, the proportion of simulations where that effect was "detected" (CI excluded zero) — this is the power estimate.

**Runtime warning:** this is an extremely expensive loop — each of the 1,000 iterations does a full `posterior_predict()` plus a complete Bayesian MCMC refit. Realistically this section alone could take many hours (or longer, depending on data size and model complexity) to finish, on top of the main model's own fit time. The inline comment `nsim <- 1000 # start small (test run)` is a little confusing/contradictory, since 1,000 is already described elsewhere in the same comment as the upper end of the "final" range — worth deciding deliberately whether you actually want to start smaller (e.g. 50–100) before committing to the full 1,000-iteration run, given the runtime.

## Output

Written to the same directory as the input file:

| File | Contents |
|---|---|
| `model.rds` | The fitted `brms` model object, saved via `saveRDS()` so it can be reloaded without refitting |
| `finCoordinationModel1.csv` | Row-level predicted probabilities of None/Asymmetrical/Symmetrical fin use |
| `finCoordinationModel2.csv` | Predicted probabilities aggregated to one row per fish × environment × training × flowSpeed cell |
| `finUsePowerModel1_partial.csv` | Checkpoint file, updated after each successful power-simulation iteration (useful if the loop is interrupted) |
| `finUsePowerModel1.csv` | Final power summary — one row per effect, with the estimated detection probability |

Default output directory (edit near the top of the script to change):
```
dataOutputPolyp
```

## Dependencies

R packages: `brms` (requires a working Stan backend — see brms installation docs), `dplyr`, `tidyr`, `ggplot2`, `readr`.

Loaded but not actually used by this script: `readxl`, `broom.mixed`, `simr`.

Install any missing packages with:
```r
install.packages(c("brms", "dplyr", "tidyr", "ggplot2", "readr"))
```

## Running it

No known blocking bugs, so it's safe to `source()` or `Rscript` top-to-bottom — just budget for the runtime:

1. Update the input path near the top and the `out_dir` path if your file locations differ.
2. Expect the main fit (4 chains × 4,000 iterations of a multinomial GLMM) to take a nontrivial amount of time on its own.
3. The power-analysis loop that follows is the expensive part (see the runtime warning above) — the checkpoint CSV lets you monitor progress or recover partial results if you need to stop it early.