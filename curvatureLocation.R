# ============================================================
# PELVIC CURVATURE LOCATION
# BINOMIAL GLM (CLIP-LEVEL, ONE CLIP PER FISH)
# ============================================================
#  What this script does, in plain terms:
#    Asks whether environment (Aquatic vs.Terrestrial) predicts where along the body peak dorsal
#    curvature tends to occur during a cycle -- Anterior vs.
#    Posterior. Cycle-level curvature-location bins are collapsed
#    to a clip-level binomial outcome (posterior_cycles out of
#    total_cycles), deduplicated to exactly one clip per fish, and
#    modeled with a single binomial GLM: environment as the only
#    fixed effect, no random effect.
# ============================================================
cat("\f")          # Clear the R console (form-feed character; works in RStudio).
rm(list = ls())     # Remove all objects from the current environment for a clean run.

# readr    - read_csv() for loading the input data
# dplyr    - data wrangling (mutate, filter, group_by, summarise, arrange, slice_tail, %>%)
# car      - Anova() for the Type II analysis-of-deviance table
# ggplot2  - the predicted-probability plot
library(readr)
library(dplyr)
library(car)
library(ggplot2)
# ============================================================
# LOAD LONG-FORM CURVATURE LOCATION DATA
# ============================================================
df_long <- read_csv(
  "curvatureLoc_long10.csv"
)
out_dir <- "dataOutput"
# ============================================================
# COMBINE INTO ANTERIOR VS POSTERIOR
# ============================================================
# curv_loc_bin is a finer-grained bin index for WHERE along the
# body peak curvature occurred; bins 2-3 are collapsed to
# "Anterior" and bins 4-5 to "Posterior". Any row whose bin isn't
# in {2,3,4,5} (e.g. bin 1 or bin 6) gets region = NA and is
# dropped as a "sparse edge bin" -- this script only analyzes a
# subset of the possible curvature-location bins.
df_long$region <- NA
df_long$region[
  df_long$curv_loc_bin %in% c(2,3)
] <- "Anterior"
df_long$region[
  df_long$curv_loc_bin %in% c(4,5)
] <- "Posterior"
# Remove sparse edge bins
df_long <- df_long %>%
  filter(!is.na(region))
# ============================================================
# CREATE CLIP-LEVEL SUMMARY
# ============================================================
# Collapse cycle-level rows to one row per clip: count how many
# cycles in this clip had peak curvature classified Posterior vs.
# Anterior, and compute the clip's posterior_prop = posterior /
# total. This is the standard setup for a binomial GLM with known
# trial counts (cbind(successes, failures) ~ ...) used below.
# NOTE: if a clip happened to have ALL of its cycles fall in the
# dropped "sparse edge bins" above, it would end up with
# total_cycles == 0 and posterior_prop = 0/0 = NaN here -- worth
# checking whether that occurs in your data, since a total_cycles
# of 0 is a degenerate binomial trial count. (It does not occur in
# the current curvatureLoc_long10.csv.)
df <- df_long %>%

  group_by(
    clip,
    fishID,
    environment,
    training
  ) %>%

  summarise(

    posterior_cycles = sum(region == "Posterior"),

    anterior_cycles = sum(region == "Anterior"),

    total_cycles = posterior_cycles + anterior_cycles,

    posterior_prop = posterior_cycles / total_cycles,

    .groups = "drop"
  )
# ============================================================
# DEDUPLICATE TO ONE CLIP PER FISH (KEEP LAST VIDEO)
# ============================================================
# Two fish (B05, greenwhite) have two clips apiece; every other fish has exactly one. Per
# project convention, when a fish was filmed twice only the LAST
# video should have been kept -- this step enforces that so every
# fish contributes exactly one row, matching the study design and
# making "clip" and "fish" the same unit of replication (i.e. no
# random effect for fishID is needed downstream). "Last" is
# determined by sorting each fish's clip names ascending (e.g.
# "..._000000" before "..._000001") and keeping the final one.
cat("\n==============================\n")
cat("DEDUPLICATION\n")
cat("==============================\n")
cat("Clip-level rows before dedup:", nrow(df), "\n")
df_dropped <- df %>%
  group_by(fishID) %>%
  filter(n() > 1) %>%
  arrange(clip, .by_group = TRUE) %>%
  slice_head(n = -1) %>%   # every row except the last, per fish
  ungroup()
cat("Rows dropped as superseded duplicate videos:\n")
print(as.data.frame(df_dropped[, c("clip", "fishID", "environment", "training")]))
df <- df %>%
  group_by(fishID) %>%
  arrange(clip, .by_group = TRUE) %>%
  slice_tail(n = 1) %>%
  ungroup()
cat("Clip-level rows after dedup (= one row per fish):", nrow(df), "\n")
# ============================================================
# FACTORS
# ============================================================
# Unlike the treadmill/swim/fin-use/fin-engagement scripts, this
# script does NOT relabel environment/training from 0/1 codes to
# "Aquatic"/"Terrestrial" or "Untrained"/"Trained" -- it just calls
# factor() on whatever values are already in the source CSV, which
# already stores them as text labels.
df <- df %>%
  mutate(
    environment = factor(environment),
    training    = factor(training)
  )
# ============================================================
# DATA SUMMARY
# ============================================================
cat("\n==============================\n")
cat("DATA SUMMARY\n")
cat("==============================\n")
print(summary(df$posterior_prop))
cat("\nGroup counts:\n")
print(table(df$environment, df$training))
# ============================================================
# ENVIRONMENT-ONLY BINOMIAL GLM  <-- THE RELEVANT TEST
# ============================================================
# One row per fish (see dedup step above), so no random effect is
# needed: each row is now a genuinely independent binomial trial.
mod_glm_env <- glm(

  cbind(
    posterior_cycles,
    anterior_cycles
  ) ~ environment,

  family = binomial,
  data = df
)
# ============================================================
# MODEL SUMMARY
# ============================================================
cat("\n==============================\n")
cat("ENVIRONMENT-ONLY GLM SUMMARY\n")
cat("==============================\n")
print(summary(mod_glm_env))
# ============================================================
# TYPE II ANALYSIS OF DEVIANCE
# ============================================================
cat("\n==============================\n")
cat("TYPE II ANALYSIS OF DEVIANCE\n")
cat("==============================\n")
anova_env <- Anova(
  mod_glm_env,
  type = 2
)
print(anova_env)
# ============================================================
# ODDS RATIOS + CONFIDENCE INTERVALS
# ============================================================
cat("\n==============================\n")
cat("ODDS RATIOS\n")
cat("==============================\n")
coef_env <- data.frame(

  effect = rownames(summary(mod_glm_env)$coefficients),

  estimate = summary(mod_glm_env)$coefficients[, "Estimate"],

  std_error = summary(mod_glm_env)$coefficients[, "Std. Error"],

  z_value = summary(mod_glm_env)$coefficients[, "z value"],

  p_value = summary(mod_glm_env)$coefficients[, "Pr(>|z|)"]
)
coef_env$odds_ratio <- exp(coef_env$estimate)
coef_env$lower_CI <- exp(
  coef_env$estimate - 1.96 * coef_env$std_error
)
coef_env$upper_CI <- exp(
  coef_env$estimate + 1.96 * coef_env$std_error
)
print(coef_env)
write.csv(
  coef_env,
  file.path(out_dir, "curvatureLocation_environmentModel_oddsRatios.csv"),
  row.names = FALSE
)
# ============================================================
# PREDICTED PROBABILITIES
# ============================================================
# SE --
# which is on the PROBABILITY scale -- and fed it directly into
# qlogis(prob) +/- 1.96*SE as if that SE were already on the LOGIT
# scale. Those two scales aren't interchangeable (a probability-
# scale SE around 0.08 is nowhere near the corresponding logit-
# scale SE around 0.48 for a probability near 0.36-0.62), so that
# produced a badly-too-narrow interval -- running the original
# formula on this deduplicated data gives Aquatic = 36% [33%, 40%],
# which is implausibly tight for n=24. The standard delta-method
# fix is to convert the probability-scale SE to the logit scale
# first (dividing by p(1-p), the derivative of the logit link at p)
# before applying it on the logit scale -- equivalent to what
# predict(type = "link", se.fit = TRUE) would give directly. With
# that fix, Aquatic = 36% [22%, 53%], a far more realistic interval
# for this sample size.
newdata_env <- expand.grid(
  environment = levels(df$environment)
)
pred_env <- predict(
  mod_glm_env,
  newdata = newdata_env,
  type = "response",
  se.fit = TRUE
)
newdata_env$predicted_prob <- pred_env$fit
newdata_env$SE <- pred_env$se.fit
newdata_env$SE_link <- newdata_env$SE / (newdata_env$predicted_prob * (1 - newdata_env$predicted_prob))
newdata_env$lower_CI <- plogis(
  qlogis(newdata_env$predicted_prob) - 1.96 * newdata_env$SE_link
)
newdata_env$upper_CI <- plogis(
  qlogis(newdata_env$predicted_prob) + 1.96 * newdata_env$SE_link
)
print(newdata_env)
write.csv(
  newdata_env,
  file.path(out_dir, "curvatureLocation_environmentModel_predictedProbabilities.csv"),
  row.names = FALSE
)
# ============================================================
# COLORS
# ============================================================
env_colors <- c(
  "Terrestrial" = "#8B4513",
  "Aquatic"     = "#1E4E8C"
)
# ============================================================
# PLOT
# ============================================================
# Bar chart of the environment-only GLM's predicted probability of
# posterior peak curvature, with 95% CI error bars. Unlike the
# original script's p_env (print()-only), this version IS saved to
# disk via ggsave() below.
p_env <- ggplot(
  newdata_env,
  aes(
    x = environment,
    y = predicted_prob,
    fill = environment
  )
) +

  geom_col(
    width = 0.7,
    alpha = 0.8
  ) +

  geom_errorbar(
    aes(
      ymin = lower_CI,
      ymax = upper_CI
    ),
    width = 0.2,
    linewidth = 0.8
  ) +

  scale_fill_manual(values = env_colors) +

  ylim(0,1) +

  labs(
    title = "Probability of Pelvic-region Peak Curvature",
    x = "Environment",
    y = "Predicted Probability"
  ) +

  theme_classic(base_size = 14) +

  theme(
    legend.position = "none"
  )
print(p_env)
# ============================================================
# SAVE
# ============================================================
ggsave(
  file.path(out_dir, "predictedProb_environment.svg"),
  p_env,
  width = 7,
  height = 5,
  dpi = 300
)
# plot mod_glm_env for diagnostics
# Base R's plot.glm() 4-panel diagnostic plot (residuals vs.
# fitted, QQ, scale-location, residuals vs. leverage) -- printed to
# the active graphics device only, not saved (same as the original
# script).
plot(mod_glm_env)
# ============================================================
# POWER ANALYSIS (ENVIRONMENT-ONLY GLM, DEDUPLICATED DATA)
# ============================================================
# frequentist parametric-bootstrap approach as the original script's
# power section (and the fin-engagement script): simulate new
# posterior_cycles counts per fish from the fitted model's own
# predicted probability, holding each fish's real total_cycles
# fixed, refit the same GLM on the simulated counts, and check
# whether the environment term's p-value clears .05. This tests
# detectability only at mod_glm_env's own fitted effect size (an
# odds ratio of ~2.83), not a swept range of hypothetical effect
# sizes, and now correctly reflects the deduplicated n=24 dataset
# and the corrected model, rather than the original script's
# n=26/pseudoreplication-affected version.
# -----------------------------
# SETTINGS
# -----------------------------
nsim <- 1000   # plain glm refits are fast, so 1000 should run quickly
set.seed(123)
# -----------------------------
# HELPER: pull a p-value for a given term from a glm summary
# -----------------------------
get_p_glm <- function(fit, term_name) {
  coefs <- summary(fit)$coefficients
  if (!(term_name %in% rownames(coefs))) return(NA)
  coefs[term_name, "Pr(>|z|)"]
}
# -----------------------------
# HELPER: binomial CI + tidy row for a power estimate
# -----------------------------
summarize_power <- function(sig_vec) {
  sig_vec <- sig_vec[!is.na(sig_vec)]
  n_valid <- length(sig_vec)
  n_sig   <- sum(sig_vec)
  p_hat   <- n_sig / n_valid
  ci      <- binom.test(n_sig, n_valid)$conf.int
  c(power = p_hat, lowerCI = ci[1], upperCI = ci[2], nsim = n_valid)
}
# -----------------------------
# SIMULATE + REFIT
# -----------------------------
total_cycles <- df$total_cycles
fitted_probs_env <- fitted(mod_glm_env)
sig_env_only <- rep(NA, nsim)
for (i in 1:nsim) {

  if (i %% 100 == 0) cat("[env-only model] Simulation", i, "/", nsim, "\n")

  posterior_sim <- rbinom(nrow(df), size = total_cycles, prob = fitted_probs_env)
  anterior_sim  <- total_cycles - posterior_sim

  fit_sim <- try(
    glm(
      cbind(posterior_sim, anterior_sim) ~ environment,
      family = binomial,
      data = df
    ),
    silent = TRUE
  )

  if (inherits(fit_sim, "try-error")) next

  p_env_sim <- get_p_glm(fit_sim, grep("^environment", rownames(summary(fit_sim)$coefficients), value = TRUE)[1])
  sig_env_only[i] <- !is.na(p_env_sim) && p_env_sim < 0.05
}
env_only_stats <- summarize_power(sig_env_only)
# ============================================================
# COMBINE + SAVE
# ============================================================
power_summary_curv <- data.frame(
  variable = "curv_loc_bin",
  model    = "environment_only_deduplicated",
  method   = "power_GLM",
  effect   = "environment",
  power    = env_only_stats["power"],
  lowerCI  = env_only_stats["lowerCI"],
  upperCI  = env_only_stats["upperCI"],
  nsim     = env_only_stats["nsim"]
)
write.csv(
  power_summary_curv,
  file.path(out_dir, "curvatureLocationPower_deduplicated.csv"),
  row.names = FALSE
)
print(power_summary_curv)
