# ============================================================
#  FIN-USE ANALYSIS — MODEL 1 ONLY (BAYESIAN MULTINOMIAL GLMM)
# ============================================================
#  What this script does, in plain terms:
#    A trimmed version of the full fin-use/coordination script that
#    keeps only Model 1: a Bayesian multinomial mixed model (via
#    brms) of finUse (a 3-level outcome: "None" / "Asymmetrical" /
#    "Symmetrical") as a function of environment, training, and
#    flowSpeed (plus their pairwise interactions), with a random
#    intercept per fish. Fit with brms::brm() using a categorical
#    likelihood, sampled via Stan/MCMC. This script:
#      - Fits the model and saves it to disk.
#      - Runs standard Bayesian diagnostics (summary, trace plots,
#        posterior predictive check).
#      - Computes predicted probabilities of each fin-use category,
#        both per-row and aggregated per fish/condition, and saves
#        both as CSVs.
#      - Plots predicted probability vs. flow speed, faceted by
#        environment x training.
#      - Runs a large simulation-based power analysis for every
#        fixed effect in the model.
#
#  WHAT WAS REMOVED, AND WHY:
#    The original version of this script also included a second,
#    separate model ("Model 2"): a frequentist binomial GLMM
#    (lme4::glmer) restricted to trials where fins were actually
#    used, asking whether use was Symmetrical vs. Asymmetrical, plus
#    its own predicted-probability table and two plots. That whole
#    section (and the `lme4` library, which was only ever used by
#    that model) has been removed here, since Model 1 is the only
#    part actually being used. Two side effects of that removal
#    worth knowing:
#      1. One of the two Model 2 plots had a missing "+" before a
#         `library(brms)` line, which would throw a runtime error and
#         halt a top-to-bottom run of the *original* script before
#         ever reaching the power analysis below. That bug no longer
#         exists in this file, because the code it was in is gone.
#      2. `coord_table` (Model 2's predicted-probability table) was
#         never written to a CSV in the original script either — so
#         nothing that used to reach disk from Model 2 is lost by
#         removing it.
#    If you ever want Model 2 back, the fuller two-model version
#    (with its bugs documented but not fixed) is the other file
#    already in this conversation.
#
#  REMAINING KNOWN ISSUES (unchanged from the full version — see
#  inline flags at each location):
#    1. `readxl` is loaded (and the section header below says "LOAD
#       DATA (EXCEL)") but the actual data load uses `read_csv()` —
#       likely leftover from an earlier version of this script that
#       read an .xlsx file directly.
#    2. `broom.mixed` and `simr` are loaded but never used anywhere
#       in this script (the power analysis here is fully custom,
#       not built on simr's powerSim()).
#    3. `library(brms)` is loaded again, redundantly, partway
#       through the file (inside the power-analysis section) — it
#       was already loaded at the top; harmless, just unnecessary.
#
#  NOTE: Aside from removing the Model 2 section and the `lme4`
#  import (its only user), no other analysis logic has been changed
#  from the original script. Comments were added to explain intent;
#  nothing else was reordered, renamed, recalculated, or fixed.
# ============================================================

cat("\f")          # Clear the R console (form-feed character; works in RStudio).
rm(list = ls())     # Remove all objects from the current environment for a clean run.

# ============================================================
# PACKAGES
# ============================================================
# brms        - Bayesian multilevel/regression models via Stan; used
#               here for the multinomial (categorical) GLMM (Model 1)
# dplyr       - data wrangling (mutate, filter, select, %>%, etc.)
# readxl      - Excel file reading — loaded but NOT used below (see
#               "known issue" #1 above; data is actually loaded from CSV)
# tidyr       - pivot_longer(), used to reshape predictions for plotting
# ggplot2     - all plotting in this script
# broom.mixed - tidy() for mixed models — loaded but not called
#               anywhere in this script (see "known issue" #2)
# readr       - read_csv() for loading the input data
# simr        - simulation-based power analysis for lme4 models —
#               loaded but not used; this script's power analysis is
#               a fully custom brms/posterior-based routine instead
#               (see "known issue" #2 and the POWER ANALYSIS section)
#
# (lme4 has been removed from this trimmed version — it was only
# ever used by Model 2's glmer() call, which is no longer here.)
library(brms)
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(broom.mixed)
library(readr)
library(simr)
# ============================================================
# LOAD DATA (EXCEL)
# ============================================================
# (Section header says "EXCEL", but this actually reads a CSV — see
# "known issue" #1 above.) One row per trial/timepoint, matching the
# treadmill dataset in structure (same fishID/environment/training/
# flowSpeed columns), plus a fin-use classification per row.
df <- read_csv("combinedWalkDataSend.csv")
# ============================================================
# OUTPUT DIRECTORY
# ============================================================
out_dir <- "dataOutputPolyp"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)   # create it if it doesn't already exist
# ============================================================
# CLEAN + FORMAT
# ============================================================
# Coerce/relabel columns for modeling:
#   - fishID: factor, for use as the mixed-model grouping variable
#   - fileName: a readable alias for the existing base_id column
#     (kept around so predictions can later be traced back to the
#     specific video/trial file they came from, e.g. in fish_means)
#   - environment / training: recoded from 0/1 codes to readable
#     factor labels (0 = Aquatic/Untrained, 1 = Terrestrial/Trained)
#   - finUse: 3-level outcome factor, with "None" set as the first
#     (and therefore reference) level — this matters below, since
#     brms's categorical() family builds separate coefficients for
#     every level relative to whichever level is first ("None" here,
#     hence the "muAsymmetrical" / "muSymmetrical" parameter names
#     used throughout the rest of the script)
#   - flowSpeed: cast to factor here, then immediately cast back to
#     numeric right after the pipe (see next line) — a round-trip
#     conversion (numeric -> factor -> character -> numeric) that
#     nets out to the same numeric values for normal flow-speed data,
#     but is redundant as written; likely leftover from earlier code
#     that treated flowSpeed categorically.
df <- df %>%
  mutate(
    fishID = factor(fishID),

    fileName = base_id,

    environment = factor(environment,
                         levels = c(0,1),
                         labels = c("Aquatic","Terrestrial")),

    training = factor(training,
                      levels = c(0,1),
                      labels = c("Untrained","Trained")),

    finUse = factor(finUse,
                    levels = c("None","Asymmetrical","Symmetrical")),

    flowSpeed = as.factor(flowSpeed)
  )
df$flowSpeed <- as.numeric(as.character(df$flowSpeed))
# ============================================================
# MULTINOMIAL GLMM (MAIN MODEL)
# ============================================================
# Bayesian multinomial logistic mixed model:
#   finUse ~ (environment + training + flowSpeed)^2 + (1 | fishID)
# The "^2" shorthand expands to all three main effects plus all
# three pairwise interactions (environment:training,
# environment:flowSpeed, training:flowSpeed). "(1 | fishID)" is a
# per-fish random intercept for repeated measures.
#
# family = categorical() fits a multinomial logistic model: with
# "None" as the reference level, brms estimates one set of
# coefficients for the log-odds of "Asymmetrical" vs "None" (prefixed
# muAsymmetrical_) and another for "Symmetrical" vs "None" (prefixed
# muSymmetrical_) — this is why coefficient names throughout the rest
# of the script come in matched Asym/Sym pairs.
#
# prior: weakly informative normal(0, 2) priors on all population-
# level ("b" class) coefficients for both non-reference categories.
# The inline comment "fix separation" refers to a common problem in
# categorical/logistic regression where a predictor near-perfectly
# separates outcome categories, causing unconstrained (frequentist)
# estimates to blow up toward +/-Infinity; a modest prior keeps
# estimates finite and regularized.
#
# control: adapt_delta = 0.99 and max_treedepth = 15 are standard
# Stan/brms sampler tweaks (higher than the ~0.8/10 defaults) used to
# reduce divergent transitions and treedepth warnings in a model this
# complex (multinomial + interactions + random effects).
#
# chains = 4, cores = 4, iter = 4000: 4 parallel MCMC chains, 4000
# iterations each (half warmup by brms default), for the main/
# reported model — a substantially more thorough fit than the
# faster settings used later for the per-simulation refits in the
# power analysis (2 chains, 2000 iter).
model <- brm(
  finUse ~ (environment + training + flowSpeed)^2 + (1 | fishID),
  data = df,
  family = categorical(),

  #  ADD PRIORS (fix separation)
  prior = set_prior("normal(0, 2)", class = "b", dpar = "muAsymmetrical") +
    set_prior("normal(0, 2)", class = "b", dpar = "muSymmetrical"),

  #  ADD SAMPLER CONTROLS (fix treedepth + divergence)
  control = list(adapt_delta = 0.99, max_treedepth = 15),

  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123
)
saveRDS(model,
        file = file.path(out_dir, "model.rds"))   # persist the fitted model object (Stan fit can be slow to reproduce)
summary(model)      # posterior estimates, credible intervals, and sampler diagnostics (Rhat, ESS)
plot(model, ask = FALSE)   # trace plots + posterior density plots for each parameter (chain mixing / convergence check)
# ============================================================
# MODEL CHECK (OPTIONAL)
# ============================================================
pp_check(model)   # posterior predictive check: compares the observed finUse distribution to draws simulated from the fitted model
# ============================================================
# PREDICTED PROBABILITIES (MULTINOMIAL)
# ============================================================
# fitted() with re_formula = NA computes POPULATION-LEVEL predicted
# probabilities (i.e. averaging over/ignoring the fishID random
# effect) for each row of the original data — summary = TRUE returns
# the posterior mean (plus SE/CI, though only the estimate columns
# are kept below) for P(None), P(Asymmetrical), P(Symmetrical) per row.
pred <- fitted(model, newdata = df, summary = TRUE, re_formula = NA)
pred_df <- cbind(df, as.data.frame(pred))
# ============================================================
# FINAL TABLE (THIS IS YOUR OUTPUT)
# ============================================================
# Keep just the identifying columns plus the three predicted-
# probability columns, renamed to short, CSV-friendly names, and
# rounded to 3 decimal places for readability.
final_table <- pred_df %>%
  dplyr::select(
    fishID,
    fileName,
    environment,
    training,
    flowSpeed,
    `Estimate.P(Y = None)`,
    `Estimate.P(Y = Asymmetrical)`,
    `Estimate.P(Y = Symmetrical)`
  ) %>%
  dplyr::rename(
    prob_None = `Estimate.P(Y = None)`,
    prob_Asym = `Estimate.P(Y = Asymmetrical)`,
    prob_Sym  = `Estimate.P(Y = Symmetrical)`
  )
final_table <- final_table %>%
  mutate(across(starts_with("prob_"), ~ round(.x, 3)))
# ============================================================
# SAVE CSV
# ============================================================
# Row-level predicted probabilities (one row per original trial/timepoint).
# NOTE: despite the filename, this is still Model 1's output — the
# "Model1"/"Model2" in these two CSV names refers to the two output
# TABLES from this one brms model (row-level vs. fish-aggregated),
# not to the removed glmer "Model 2".
write.csv(final_table,
          file.path(out_dir, "finCoordinationModel1.csv"),
          row.names = FALSE)
# Aggregate to one row per fish x environment x training x flowSpeed
# cell (averaging predicted probabilities across trials/files within
# that cell), keeping a semicolon-joined list of the contributing
# file names for traceability back to the raw trials.
fish_means <- final_table %>%
  group_by(fishID, environment, training, flowSpeed) %>%
  summarise(
    fileNames = paste(unique(fileName), collapse = "; "),
    prob_None = mean(prob_None),
    prob_Asym = mean(prob_Asym),
    prob_Sym  = mean(prob_Sym),
    .groups = "drop"
  )
write.csv(fish_means,
          file.path(out_dir, "finCoordinationModel2.csv"),
          row.names = FALSE)
# ============================================================
# PLOT MULTINOMIAL
# ============================================================
# Reshape final_table from wide (one column per category) to long
# (one row per category per original row) so all three probability
# curves can share one ggplot color aesthetic, then plot predicted
# probability vs flowSpeed, faceted by environment x training, one
# line per fin-use category.
fin_colors <- c(
  "Symmetrical"  = "#3B9AB2",  # teal
  "Asymmetrical" = "#F21A00",  # orange
  "None"         = "#D0D0D0" # light grey
)
plot_df <- final_table %>%
  pivot_longer(
    cols = c(prob_None, prob_Asym, prob_Sym),
    names_to = "finUse",
    values_to = "prob"
  )
plot_df$flowSpeed <- as.numeric(as.character(plot_df$flowSpeed))
plot_df$finUse <- recode(plot_df$finUse,
                         prob_None = "None",
                         prob_Asym = "Asymmetrical",
                         prob_Sym  = "Symmetrical")
ggplot(plot_df,
       aes(x = flowSpeed, y = prob,
           color = finUse,
           group = finUse)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  facet_grid(environment ~ training) +
  scale_color_manual(values = fin_colors) +
  labs(
    x = "Flow speed (BL/s)",
    y = "Probability",
    color = "Fin use"
  ) +
  theme_classic(base_size = 14)
# ============================================================
# POWER ANALYSIS (MULTINOMIAL MODEL)
# ============================================================
# Everything from here down analyzes Model 1 (the brms multinomial
# `model` object) — this is the only model in this trimmed script.
#
# Unlike the treadmill/swim scripts' power analyses — which inject a
# range of hypothetical effect sizes (1x/2x/3x an LMM estimate) to
# ask "how big would the effect need to be to reliably detect it?" —
# this is a *retrospective / self-consistency* power check: it asks
# "if the world really looks like what this fitted model says
# (i.e., these exact posterior estimates are the truth), and I
# collected a new dataset of the same size and refit, how often would
# each effect's 95% credible interval exclude zero?" It tests
# detectability at the model's own observed effect size only, not at
# multiple hypothetical sizes.
library(brms)
# -----------------------------
# FUNCTION: check if main-effect / interaction coefficient detected
# -----------------------------
# "Detected" here means: the 95% credible interval for this
# coefficient (2.5th-97.5th percentile of its posterior draws)
# excludes zero. Returns NA if the named parameter isn't found in
# the model's posterior draws (e.g. a term dropped for some reason).
check_effect <- function(fit, param_name) {

  # extract posterior draws directly
  draws <- as_draws_df(fit)

  # check parameter exists
  if (!(param_name %in% colnames(draws))) {
    return(NA)
  }

  vals <- draws[[param_name]]

  ci <- quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE)

  return(as.numeric(ci[1] > 0 | ci[2] < 0))
}
# -----------------------------
# FUNCTION: check if a simple slope (main + interaction) is detected
# -----------------------------
# For a categorical predictor's reference level (Aquatic, here), the
# flowSpeed main-effect coefficient IS the flowSpeed slope. For the
# non-reference level (Terrestrial), the effective slope is the main
# effect PLUS the environment:flowSpeed interaction coefficient — this
# function sums those two posterior draw vectors elementwise (so the
# combination is evaluated posterior-draw-by-posterior-draw, properly
# propagating uncertainty) and checks whether that combined effect's
# 95% CI excludes zero, i.e. whether flowSpeed has a detectable slope
# specifically WITHIN terrestrial-reared fish.
check_simple_slope <- function(fit, main, interaction) {
  draws <- as_draws_df(fit)
  if (!all(c(main, interaction) %in% colnames(draws))) {
    return(NA)
  }
  combined <- draws[[main]] + draws[[interaction]]
  ci <- quantile(combined, probs = c(0.025, 0.975), na.rm = TRUE)
  return(as.numeric(ci[1] > 0 | ci[2] < 0))
}
# -----------------------------
# SETTINGS
# -----------------------------

nsim <- 1000   #  start small (test run)
# increase to 500-1000 for final
# Pre-allocate one row per planned simulation, with a column for
# every main effect / interaction / simple-slope being power-tested
# (all initialized to NA; only successful simulations fill these in —
# see the loop below).
results <- data.frame(
  sim = 1:nsim,

  # main effects
  env_asym = NA, env_sym = NA,
  train_asym = NA, train_sym = NA,
  flow_asym = NA, flow_sym = NA,

  # interactions
  env_train_asym = NA, env_train_sym = NA,
  env_flow_asym = NA, env_flow_sym = NA,
  train_flow_asym = NA, train_flow_sym = NA,

  # terrestrial simple slope (flow speed, within terrestrial fish)
  terr_slope_asym = NA, terr_slope_sym = NA
)
df_model <- model$data   # the exact data frame brms used to fit Model 1 (same rows/columns df was reduced to internally)
set.seed(123)
# -----------------------------
# SIMULATION LOOP
# -----------------------------
# For each of nsim (1000) iterations:
#   1. Draw one full posterior-predictive replicate of finUse for
#      every row (posterior_predict() returns a [n_posterior_draws x
#      n_rows] matrix of simulated category codes; one whole draw
#      -- i.e. one row of that matrix -- is picked at random).
#   2. Replace the real finUse column with that simulated one,
#      keeping every other column (predictors, fishID) identical.
#   3. Refit the multinomial model on this simulated dataset (reusing
#      the already-compiled Stan model via recompile = FALSE for
#      speed, but with lighter sampler settings than the main fit:
#      2 chains x 2000 iterations instead of 4 x 4000, and a fresh
#      seed per iteration).
#   4. For each fixed effect (and the terrestrial flowSpeed simple
#      slope), check whether its 95% credible interval excludes zero
#      in this refit — i.e., would this effect have been "detected"
#      in a fresh sample drawn from a world where the original
#      fitted model is exactly correct?
#   5. Write a checkpoint CSV after every successful iteration, so
#      long-running progress isn't lost if the loop is interrupted.
#      NOTE: because the checkpoint write is the LAST line in the
#      loop body and a failed refit (`try-error`) hits `next` before
#      reaching it, a failed iteration does NOT trigger its own
#      checkpoint write — the file simply reflects however many
#      simulations most recently succeeded.
#
# This is an extremely expensive loop: each of the 1000 iterations
# performs a full posterior_predict() plus a complete Bayesian MCMC
# refit (2 chains x 2000 iterations of Stan sampling), so this
# section alone can realistically take many hours (or longer,
# depending on data size/model complexity) to finish.
for (i in 1:nsim) {

  cat("Simulation", i, "\n")

  # ---- simulate response ----
  sim_y <- posterior_predict(model)

  # randomly select ONE simulated dataset
  draw_id <- sample(1:nrow(sim_y), 1)
  sim_vec <- sim_y[draw_id, ]

  df_sim <- df_model

  df_sim$finUse <- factor(
    levels(df_model$finUse)[sim_vec],
    levels = levels(df_model$finUse)
  )

  if (i == 1) {
    print(table(df_sim$finUse))
  }

  # ---- refit model (reuse compiled model for speed) ----
  fit_sim <- try(
    update(
      model,
      newdata = df_sim,
      recompile = FALSE,
      chains = 2,
      iter = 2000,
      cores = 4,
      refresh = 0,
      seed = 123 + i
    ),
    silent = TRUE
  )

  if (inherits(fit_sim, "try-error")) {
    cat("Model failed at simulation", i, "\n")
    next
  }

  # -----------------------------
  # MAIN EFFECTS
  # -----------------------------
  results$env_asym[i] <- check_effect(fit_sim, "b_muAsymmetrical_environmentTerrestrial")
  results$env_sym[i]  <- check_effect(fit_sim, "b_muSymmetrical_environmentTerrestrial")

  results$train_asym[i] <- check_effect(fit_sim, "b_muAsymmetrical_trainingTrained")
  results$train_sym[i]  <- check_effect(fit_sim, "b_muSymmetrical_trainingTrained")

  results$flow_asym[i] <- check_effect(fit_sim, "b_muAsymmetrical_flowSpeed")
  results$flow_sym[i]  <- check_effect(fit_sim, "b_muSymmetrical_flowSpeed")

  # -----------------------------
  # INTERACTIONS
  # -----------------------------
  results$env_train_asym[i] <- check_effect(fit_sim, "b_muAsymmetrical_environmentTerrestrial:trainingTrained")
  results$env_train_sym[i]  <- check_effect(fit_sim, "b_muSymmetrical_environmentTerrestrial:trainingTrained")

  results$env_flow_asym[i] <- check_effect(fit_sim, "b_muAsymmetrical_environmentTerrestrial:flowSpeed")
  results$env_flow_sym[i]  <- check_effect(fit_sim, "b_muSymmetrical_environmentTerrestrial:flowSpeed")

  results$train_flow_asym[i] <- check_effect(fit_sim, "b_muAsymmetrical_trainingTrained:flowSpeed")
  results$train_flow_sym[i]  <- check_effect(fit_sim, "b_muSymmetrical_trainingTrained:flowSpeed")

  # -----------------------------
  # TERRESTRIAL SIMPLE SLOPE (flow speed, within terrestrial fish)
  # -----------------------------
  results$terr_slope_asym[i] <- check_simple_slope(
    fit_sim,
    "b_muAsymmetrical_flowSpeed",
    "b_muAsymmetrical_environmentTerrestrial:flowSpeed"
  )
  results$terr_slope_sym[i] <- check_simple_slope(
    fit_sim,
    "b_muSymmetrical_flowSpeed",
    "b_muSymmetrical_environmentTerrestrial:flowSpeed"
  )

  # ---- checkpoint: write partial results after every simulation ----
  write.csv(results, file.path(out_dir, "finUsePowerModel1_partial.csv"), row.names = FALSE)
}
# -----------------------------
# POWER SUMMARY (FINAL OUTPUT)
# -----------------------------
# Collapse the per-simulation 0/1/NA detection results into one power
# estimate per effect: the proportion of (non-NA) simulations in
# which that effect's 95% credible interval excluded zero.
power_summary <- data.frame(
  effect = c(
    "Environment (Asym)", "Environment (Sym)",
    "Training (Asym)", "Training (Sym)",
    "FlowSpeed (Asym)", "FlowSpeed (Sym)",
    "Env×Train (Asym)", "Env×Train (Sym)",
    "Env×Flow (Asym)", "Env×Flow (Sym)",
    "Train×Flow (Asym)", "Train×Flow (Sym)",
    "Terrestrial simple slope, flow speed (Asym)",
    "Terrestrial simple slope, flow speed (Sym)"
  ),

  power = c(
    mean(results$env_asym, na.rm = TRUE),
    mean(results$env_sym, na.rm = TRUE),

    mean(results$train_asym, na.rm = TRUE),
    mean(results$train_sym, na.rm = TRUE),

    mean(results$flow_asym, na.rm = TRUE),
    mean(results$flow_sym, na.rm = TRUE),

    mean(results$env_train_asym, na.rm = TRUE),
    mean(results$env_train_sym, na.rm = TRUE),

    mean(results$env_flow_asym, na.rm = TRUE),
    mean(results$env_flow_sym, na.rm = TRUE),

    mean(results$train_flow_asym, na.rm = TRUE),
    mean(results$train_flow_sym, na.rm = TRUE),

    mean(results$terr_slope_asym, na.rm = TRUE),
    mean(results$terr_slope_sym, na.rm = TRUE)
  )
)
write.csv(power_summary,
          file.path(out_dir, "finUsePowerModel1.csv"),
          row.names = FALSE)
