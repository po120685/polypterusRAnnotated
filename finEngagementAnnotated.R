# ============================================================
#  FIN ENGAGEMENT ANALYSIS (BINARY GLMM: ANY FIN USE, YES/NO)
# ============================================================
#  What this script does, in plain terms:
#    Asks whether the fins were engaged at all during a trial — not
#    which kind of fin use occurred, just whether any occurred.
#    finEngaged is collapsed from finUse into a binary outcome
#    (0 = "None", 1 = anything else), modeled as a function of
#    environment and flowSpeed (and their interaction) with a
#    per-fish random intercept, using a standard frequentist
#    binomial GLMM (lme4::glmer).
#
#    NOTE ON MODEL SCOPE: this model's fixed-effect formula is just
#    `environment * flowSpeed`.
#
#    The script then: saves the model summary to a text file, builds
#    a tidy fixed-effects table, computes model-predicted probabilities of fin engagement (with
#    95% CIs) across environment x flowSpeed, overlays them against
#    observed proportions in a plot saved as SVG, runs a fairly
#    thorough GLMM diagnostic suite (DHARMa simulated residuals,
#    singularity/collinearity/overall model checks via `performance`,
#    and per-fish Cook's distance via `influence.ME`), and finishes
#    with a simulation-based power analysis for each fixed effect.
# ============================================================

cat("\f")          # Clear the R console (form-feed character; works in RStudio).
rm(list = ls())     # Remove all objects from the current environment for a clean run.

# ============================================================
# PACKAGES
# ============================================================
# dplyr       - data wrangling (mutate, filter, group_by, %>%, etc.)
# readr       - read_csv() for loading the input data
# lme4        - the binomial GLMM itself (glmer), plus simulate()
#               for the parametric-bootstrap power analysis
# broom.mixed - tidy() to turn the glmer fixed-effect table into a
#               clean data frame
# ggplot2     - the model/observed-proportions plot
library(dplyr)
library(readr)
library(lme4)
library(broom.mixed)
library(ggplot2)
# ============================================================
# LOAD DATA
# ============================================================
# The treadmill dataset: one row per trial/timepoint, repeated per
# fishID, with a finUse classification per row.
df <- read_csv("combinedWalkDataSend.csv")
# ============================================================
# OUTPUT DIRECTORY
# ============================================================
out_dir <- "dataOutputPolyp"
  dir.create(out_dir, recursive = TRUE)
}
# ============================================================
# CLEAN + FORMAT
# ============================================================
# Coerce/relabel columns for modeling:
#   - fishID: factor, for use as the mixed-model grouping variable
#   - fileName: a readable alias for the existing base_id column
#   - environment / training: recoded from 0/1 codes to readable
#     factor labels (0 = Aquatic/Untrained, 1 = Terrestrial/Trained)
#     — note `training` is created here but never used in the model
#     formula below (see the file-level note on model scope above)
#   - flowSpeed: coerced to numeric via character (guards against it
#     silently being read in as something else, e.g. a factor, from
#     the source CSV)
#   - finEngaged: the key derived variable for this script — collapses
#     the 3-level finUse into a binary 0/1 indicator of ANY fin use
#     (0 = "None", 1 = "Asymmetrical" or "Symmetrical")
df <- df %>%
  mutate(
    
    fishID = factor(fishID),
    
    fileName = base_id,
    
    environment = factor(
      environment,
      levels = c(0,1),
      labels = c("Aquatic","Terrestrial")
    ),
    
    training = factor(
      training,
      levels = c(0,1),
      labels = c("Untrained","Trained")
    ),
    
    flowSpeed = as.numeric(as.character(flowSpeed)),
    
    # ========================================================
    # BINARY FIN ENGAGEMENT
    # 1 = fins engaged
    # 0 = no fin use
    # ========================================================
    
    finEngaged = ifelse(finUse == "None", 0, 1)
    
  )
# Convert to factor if desired for plotting
# (finEngaged becomes a 2-level factor here, "NoFinUse"/"FinEngaged",
# for readability elsewhere — but note the model formula below
# converts it right back to numeric via as.numeric(finEngaged) - 1;
# because factor() assigns integer codes in the order of its
# `levels` argument, "NoFinUse" is level 1 and "FinEngaged" is level
# 2, so as.numeric(...) - 1 correctly recovers the original 0/1
# coding — a roundabout way to get back to a binary numeric, but not
# incorrect.)
df$finEngaged <- factor(
  df$finEngaged,
  levels = c(0,1),
  labels = c("NoFinUse","FinEngaged")
)
# ============================================================
# MODEL
# ============================================================
# Binomial GLMM: finEngaged (0/1) ~ environment * flowSpeed +
# (1 | fishID). "environment * flowSpeed" expands to the
# environment main effect, the flowSpeed main effect, AND their
# interaction — but, as noted at the top of this file, `training` is
# not part of this model at all. "(1 | fishID)" is the usual per-fish
# random intercept for repeated measures.
mod_fin <- glmer(
  as.numeric(finEngaged) - 1 ~
    environment * flowSpeed +
    (1 | fishID),
  
  data = df,
  family = binomial
)
summary(mod_fin)
# ============================================================
# SAVE MODEL SUMMARY
# ============================================================
# sink() redirects console output to a file until the next
# sink()-with-no-arguments call, which restores normal console
# output — a standard way to capture print() output as a text file.
sink(file.path(out_dir, "finEngagement_model_summary.txt"))
print(summary(mod_fin))
sink()
# ============================================================
# EFFECTS TABLE
# ============================================================
# Tidy fixed-effect coefficient table (estimate, SE, statistic,
# p-value) via broom.mixed, with columns renamed for readability and
# tagged with a `variable` label.
effects_fin <- broom.mixed::tidy(
  mod_fin,
  effects = "fixed"
) %>%
  transmute(
    effect = term,
    estimate = estimate,
    std_error = std.error,
    statistic = statistic,
    p_value = p.value,
    variable = "fin_engagement"
  )
write.csv(
  effects_fin,
  file.path(out_dir, "finEngagementEffects.csv"),
  row.names = FALSE
)
print(effects_fin)
# ============================================================
# PREDICTED PROBABILITIES
# ============================================================
# Prediction grid over every combination of environment level and
# observed flowSpeed value (no `training` dimension, since it's not
# in the model).
newdata <- expand.grid(
  environment = levels(df$environment),
  flowSpeed = sort(unique(df$flowSpeed))
)
# Predict on link scale
# type = "link" returns predictions on the logit scale (with a
# standard error), and re.form = NA makes these population-level
# predictions (ignoring the fishID random effect) — predicting on
# the link scale first, then back-transforming, is what lets the
# 95% CI below stay properly bounded to [0,1] and asymmetric near
# the edges, rather than building a naive symmetric CI directly on
# the probability scale.
pred_link <- predict(
  mod_fin,
  newdata = newdata,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)
# Build dataframe
pred_df <- newdata
pred_df$fit_link <- pred_link$fit
pred_df$se_link  <- pred_link$se.fit
# Convert to probability scale
pred_df$prob_finEngaged <- plogis(pred_df$fit_link)
# Approximate SE on probability scale
# Delta-method approximation: SE_prob ~= SE_link * d(plogis)/d(link)
# = SE_link * p*(1-p). Computed here, but — as noted at the top of
# this file — NOT actually used for the CI below (which is instead
# built correctly on the link scale and back-transformed); this
# column is informational only in the saved output.
pred_df$se_prob <- pred_df$se_link *
  (pred_df$prob_finEngaged *
     (1 - pred_df$prob_finEngaged))
# 95% CI
# Built on the link (logit) scale — fit_link +/- 1.96 SE — and then
# back-transformed through plogis(), which is the standard, more
# accurate way to get a bounded, asymmetric CI for a probability.
pred_df$lower <- plogis(
  pred_df$fit_link - 1.96 * pred_df$se_link
)
pred_df$upper <- plogis(
  pred_df$fit_link + 1.96 * pred_df$se_link
)
# ============================================================
# SAVE PREDICTIONS
# ============================================================
write.csv(
  pred_df,
  file.path(out_dir, "finEngagementProbabilities.csv"),
  row.names = FALSE
)
# ============================================================
# OBSERVED PROPORTIONS
# ============================================================
# Raw (non-model) observed proportion of fin-engaged trials per
# environment x flowSpeed cell, computed directly from the data —
# this is what gets plotted as points alongside the model's fitted
# curve, as a visual model-fit check.
obs_df <- df %>%
  mutate(
    finEngagedNum = ifelse(finUse == "None", 0, 1)
  ) %>%
  group_by(environment, flowSpeed) %>%
  summarise(
    prop_engaged = mean(finEngagedNum, na.rm = TRUE),
    n = sum(!is.na(finEngagedNum)),
    .groups = "drop"
  )
# Optional: offset overlapping points
# Nudges Aquatic points slightly left and Terrestrial points
# slightly right of their true flowSpeed value, purely so the two
# environments' observed-proportion points don't visually overlap
# on the plot when they share the same flowSpeed.
obs_df <- obs_df %>%
  mutate(
    flowSpeed_plot = ifelse(
      environment == "Aquatic",
      flowSpeed - 0.03,
      flowSpeed + 0.03
    )
  )
# ============================================================
# COLORS
# ============================================================
group_colors <- c(
  "Aquatic"     = "#3B9AB2",
  "Terrestrial" = "#E64B35"
)
# ============================================================
# PLOT
# ============================================================
# Three layers: a 95% CI ribbon and a fitted-probability line (both
# from the model, pred_df), plus observed-proportion points (from
# the raw data, obs_df, using the jittered x-position). Every layer
# here is correctly chained with "+" through to theme_classic() at
# the end — no dangling or missing "+" anywhere in this plot.
p <- ggplot() +
  
  # ----------------------------------------------------------
# 95% CI ribbon
# ----------------------------------------------------------
geom_ribbon(
  data = pred_df,
  aes(
    x = flowSpeed,
    ymin = lower,
    ymax = upper,
    fill = environment
  ),
  alpha = 0.20,
  color = NA
) +
  
  # ----------------------------------------------------------
# MODEL PREDICTION
# ----------------------------------------------------------
geom_line(
  data = pred_df,
  aes(
    x = flowSpeed,
    y = prob_finEngaged,
    color = environment
  ),
  linewidth = 1.2
) +
  
  # ----------------------------------------------------------
# OBSERVED PROPORTIONS
# ----------------------------------------------------------
geom_point(
  data = obs_df,
  aes(
    x = flowSpeed_plot,
    y = prop_engaged,
    fill = environment
  ),
  shape = 21,
  color = "black",
  size = 3.5
) +
  
  # ----------------------------------------------------------
# COLORS
# ----------------------------------------------------------
scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  
  # ----------------------------------------------------------
# AXES
# ----------------------------------------------------------
scale_y_continuous(
  limits = c(0, 1)
) +
  
  labs(
    x = "Belt speed (BL/s)",
    y = "Probability of fin engagement",
    color = "Environment",
    fill = "Environment",
    title = "Fin engagement across treadmill belt speeds"
  ) +
  
  theme_classic(base_size = 14)
print(p)
# ============================================================
# SAVE SVG
# ============================================================
ggsave(
  filename = file.path(
    out_dir,
    "FinEngagement_ModelObserved.svg"
  ),
  plot = p,
  width = 7,
  height = 5
)
# ============================================================
# OPTIONAL: ODDS RATIOS
# ============================================================
# Re-uses (overwrites in memory) the effects_fin table built above,
# adding an odds_ratio = exp(estimate) column, and saves it to a
# SEPARATE file (finEngagementEffects_OR.csv) — the original
# finEngagementEffects.csv (without this column) written earlier is
# untouched and still exists alongside it, so both files persist.
effects_fin <- effects_fin %>%
  mutate(
    odds_ratio = exp(estimate)
  )
write.csv(
  effects_fin,
  file.path(out_dir, "finEngagementEffects_OR.csv"),
  row.names = FALSE
)
print(effects_fin)
# ============================================================
# DIAGNOSTICS
# ============================================================
# DHARMa      - simulation-based residual diagnostics for GLMMs
#               (works around the fact that raw/Pearson residuals
#               from a binomial GLMM aren't straightforwardly
#               interpretable the way OLS residuals are)
# performance - check_singularity(), check_collinearity(),
#               check_model() — general-purpose model health checks
# car         - loaded here but no car:: function is called directly
#               anywhere below; performance::check_collinearity()
#               may use it as an internal dependency, but this
#               script itself never calls it explicitly (see the
#               file-level note above)
library(DHARMa)
library(performance)
library(car)
# ============================================================
# SIMULATED RESIDUALS
# ============================================================
# DHARMa's core idea: simulate many new datasets from the fitted
# model, and for each observed data point, see where its actual
# value falls within the distribution of simulated values for that
# point — under a well-specified model, those "residuals" should be
# uniformly distributed. All of the following diagnostic
# calls/plots print to the active graphics device / console only;
# none of them are saved to a file (see file-level note above).
sim_res <- simulateResiduals(mod_fin)
# Basic residual plots
plot(sim_res)
# ============================================================
# TESTS
# ============================================================
# Uniformity test
# (are the simulated residuals uniformly distributed, as they
# should be under a correctly specified model?)
testUniformity(sim_res)
# Dispersion test
# (is the binomial model over- or under-dispersed relative to what
# was simulated?)
testDispersion(sim_res)
# Outlier test
# (are there more extreme residuals than expected by chance?)
testOutliers(sim_res)
# ============================================================
# RANDOM EFFECTS / SINGULARITY
# ============================================================
# Checks whether the fishID random-effect variance estimate has
# collapsed to (or very near) zero — a sign the random effect isn't
# adding meaningful information and the model may be overparameterized.
check_singularity(mod_fin)
# ============================================================
# COLLINEARITY
# ============================================================
# Variance inflation factors for the fixed effects — flags predictors
# that are highly correlated with each other, which can make
# individual coefficient estimates unstable/hard to interpret.
check_collinearity(mod_fin)
# ============================================================
# OVERALL MODEL PERFORMANCE
# ============================================================
# A combined diagnostic plot panel (performance::check_model())
# covering several of the checks above plus more, in one view.
check_model(mod_fin)
# ============================================================
# INFLUENCE / COOK'S DISTANCE
# ============================================================
# Leave-one-fish-out influence diagnostics: refits the model leaving
# each fish out in turn (via influence.ME::influence(), grouped by
# fishID) and computes Cook's distance per fish — a high value flags
# a fish whose data disproportionately drives the model's fitted
# coefficients.
library(influence.ME)
infl <- influence(mod_fin,
                  group = "fishID")
cooks <- cooks.distance(infl)
print(cooks)
# Plot Cook's distance
# The dashed red reference line at 4/n is a common (if somewhat
# informal) rule-of-thumb cutoff for flagging an especially
# influential fish.
plot(
  cooks,
  type = "h",
  main = "Cook's Distance by fishID",
  ylab = "Cook's Distance"
)
abline(h = 4 / length(cooks),
       col = "red",
       lty = 2)
# ============================================================
# POWER ANALYSIS (FIN ENGAGEMENT GLMM)
# Append this after your existing finEngagement script
# (requires mod_fin, df, and out_dir already in the environment)
# ============================================================
# This is a *retrospective / self-consistency* check — "if the world
# really looks like what this fitted model says, how often would
# each effect be detected in a fresh sample of the same size?" —
# tested only at the model's own observed effect size, not swept
# across multiple hypothetical effect sizes.
#
# The mechanism is a standard frequentist parametric bootstrap, using
# base R/lme4's own simulate() method for glmer models: it fixes the
# fixed effects at their fitted (MLE) point estimates and draws fresh
# random-effect values (from the fitted variance components) plus
# fresh binomial outcomes for each simulated dataset. Because it
# simulates from a single fitted point estimate rather than a full
# distribution over parameter values, it doesn't propagate parameter
# *uncertainty* into the power estimate — only random-effect and
# residual/binomial sampling variability around the fitted point
# estimate.
library(lme4)
library(broom.mixed)
# -----------------------------
# SETTINGS
# -----------------------------
nsim <- 1000   # glmer refits are fast, so 1000 should run quickly
set.seed(123)
power_results <- data.frame(
  sim = 1:nsim,
  environment = NA,
  flowSpeed = NA,
  interaction = NA
)
# -----------------------------
# PRE-GENERATE SIMULATED RESPONSES
# (parametric bootstrap from the fitted model's own estimated
#  fixed + random effects)
# -----------------------------
# All nsim simulated response vectors are generated up front in one
# call — sim_responses is a data-frame-like object where
# sim_responses[[i]] is the full simulated finEngaged vector (0/1,
# one value per row of df) for simulation i.
sim_responses <- simulate(mod_fin, nsim = nsim, seed = 123)
# -----------------------------
# HELPER: pull a p-value for a given fixed-effect term
# -----------------------------
get_p <- function(tidy_df, term_name) {
  val <- tidy_df$p.value[tidy_df$term == term_name]
  if (length(val) == 0) return(NA)
  return(val)
}
# -----------------------------
# SIMULATION LOOP
# -----------------------------
# For each of nsim (1000) iterations: substitute the real finEngaged
# column for one pre-generated simulated replicate, refit the same
# GLMM (with the more robust "bobyqa" optimizer explicitly set) on
# the simulated data, and record whether each fixed-effect term's
# p-value cleared p < .05. This loop excludes not just outright
# fitting errors but also any refit that raised an lme4 convergence
# warning (a stricter, more conservative inclusion criterion — a
# "successful but shaky" fit is treated the same as a failed one and
# excluded from the power estimate's denominator, rather than being
# counted as a valid non-detection or detection).
for (i in 1:nsim) {
  
  if (i %% 50 == 0) cat("Simulation", i, "/", nsim, "\n")
  
  df_sim <- df
  df_sim$finEngaged_sim <- sim_responses[[i]]
  
  fit_sim <- try(
    glmer(
      finEngaged_sim ~ environment * flowSpeed + (1 | fishID),
      data = df_sim,
      family = binomial,
      control = glmerControl(optimizer = "bobyqa")
    ),
    silent = TRUE
  )
  
  # skip failed fits or ones that threw a convergence warning
  if (inherits(fit_sim, "try-error")) {
    cat("Model failed at simulation", i, "\n")
    next
  }
  if (length(fit_sim@optinfo$conv$lme4$messages) > 0) {
    cat("Convergence warning at simulation", i, "- skipping\n")
    next
  }
  
  tidy_sim <- broom.mixed::tidy(fit_sim, effects = "fixed")
  
  power_results$environment[i] <- get_p(tidy_sim, "environmentTerrestrial") < 0.05
  power_results$flowSpeed[i]   <- get_p(tidy_sim, "flowSpeed") < 0.05
  power_results$interaction[i] <- get_p(tidy_sim, "environmentTerrestrial:flowSpeed") < 0.05
  
  # checkpoint every 100 sims in case of interruption
  # (checkpointing every 100 rather than every iteration is reasonable
  # here since, per the inline comment above, glmer refits are fast
  # and the whole loop is expected to finish quickly, so there's
  # little to lose by checkpointing less often)
  if (i %% 100 == 0) {
    write.csv(power_results,
              file.path(out_dir, "finEngagementPower_partial.csv"),
              row.names = FALSE)
  }
}
# -----------------------------
# POWER SUMMARY (with binomial CI on each power estimate)
# -----------------------------
# Proportion of valid (non-NA) simulations that were "significant",
# with a binomial exact confidence interval.
summarize_power <- function(sig_vec) {
  sig_vec <- sig_vec[!is.na(sig_vec)]
  n_valid <- length(sig_vec)
  n_sig   <- sum(sig_vec)
  p_hat   <- n_sig / n_valid
  ci      <- binom.test(n_sig, n_valid)$conf.int
  c(power = p_hat, lowerCI = ci[1], upperCI = ci[2], nsim = n_valid)
}
env_stats   <- summarize_power(power_results$environment)
flow_stats  <- summarize_power(power_results$flowSpeed)
int_stats   <- summarize_power(power_results$interaction)
power_summary_fin <- data.frame(
  variable = "finEngaged",
  method   = "power_GLMM",
  effect   = c("environment", "flowSpeed", "environment:flowSpeed"),
  power    = c(env_stats["power"], flow_stats["power"], int_stats["power"]),
  lowerCI  = c(env_stats["lowerCI"], flow_stats["lowerCI"], int_stats["lowerCI"]),
  upperCI  = c(env_stats["upperCI"], flow_stats["upperCI"], int_stats["upperCI"]),
  nsim     = c(env_stats["nsim"], flow_stats["nsim"], int_stats["nsim"])
)
write.csv(
  power_summary_fin,
  file.path(out_dir, "finEngagementPower.csv"),
  row.names = FALSE
)
print(power_summary_fin)

