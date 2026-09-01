# ============================================================
#  TREADMILL MASTER ANALYSIS
# ============================================================
#  What this script does, in plain terms:
#    This is a heavier sibling of the terrestrial-gait ANOVA script:
#    instead of a simple two-way ANOVA per variable, it fits mixed
#    models (subjects have repeated measurements, so `fishID` is
#    included as a random intercept) and, for two of the eight
#    outcome variables, falls back to a permutation-based ANOVA
#    instead of a parametric one.
#
#    For each of 8 treadmill/kinematic variables, per a manually
#    chosen method ("LMM", "logLMM", or "perm"):
#      - LMM     : linear mixed model, y ~ (environment + training +
#                  flowSpeed)^2 + (1 | fishID), fit with lmer().
#      - logLMM  : same, but on log(y) (for right-skewed variables).
#      - perm    : a permutation ANOVA (aovperm) on the same fixed
#                  effects, clustering permutations by fishID — used
#                  when a variable doesn't fit the LMM's parametric
#                  assumptions well.
#
#    For every variable, the script also:
#      - Runs post-hoc pairwise/trend comparisons for any of the
#        three two-way interactions that came out significant.
#      - Records model diagnostics (normality of residuals,
#        heteroscedasticity) for LMM/logLMM variables.
#      - Runs a power analysis for every fixed-effect term:
#          * LMM/logLMM  -> simulation-based power via simr::powerSim()
#            at the model's own observed effect size.
#          * perm        -> a custom simulation: inject a known effect
#            of 1x/2x/3x the LMM-estimated coefficient into simulated
#            data, rerun the permutation ANOVA 1000 times per effect
#            size, and report the fraction of runs where p < .05
#            (with a binomial exact confidence interval).
#
#    Finally, results from all variables are stacked into four
#    tables (model results, diagnostics, post-hoc comparisons, power)
#    and written out as CSVs.
#
#  NOTE: This is an annotation pass only -- no analysis logic has
#  been changed from the original script. Comments were added to
#  explain intent; nothing was reordered, renamed, or recalculated.
# ============================================================

cat("\f")          # Clear the R console (form-feed character; works in RStudio).
rm(list = ls())     # Remove all objects from the current environment for a clean run.

# ============================================================
#   TREADMILL MASTER ANALYSIS
# ============================================================
# readr       - fast CSV I/O (read_csv)
# dplyr       - data wrangling (mutate, filter, %>%, etc.)
# lme4        - linear mixed-effects models (lmer, fixef)
# lmerTest    - adds p-values / Satterthwaite df to lme4 models so
#               anova(model) and summary(model) report significance
# emmeans     - estimated marginal means, pairwise contrasts
#               (pairs()), and slope comparisons (emtrends())
# performance - model diagnostics: check_normality(),
#               check_heteroscedasticity(), r2()
# permuco     - permutation-based ANOVA (aovperm), used as a
#               distribution-free alternative to the LMM
# simr        - simulation-based power analysis for mixed models
#               (powerSim, fixed())
# tidyr       - drop_na() for listwise-deleting missing predictors
library(readr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(permuco)
library(simr)
library(tidyr)

# Fixed seed so permutation ANOVAs and power simulations (both of
# which rely on random number generation) are reproducible run to run.
set.seed(123)

# ============================================================
# RESULTS STORAGE
# ============================================================
# Model coefficient/ANOVA-style results (rbind'd across every
# variable and method) end up here.
results_table <- data.frame()

# ============================================================
# DIAGNOSTICS STORAGE
# ============================================================
# Pre-declared with explicit column types/names so rbind() behaves
# predictably even before any rows are added (LMM diagnostics get
# real values; perm-method rows get NA — see run_perm()).
diagnostics_table <- data.frame(
  variable = character(),
  method = character(),
  normality_test = character(),
  normality_p = numeric(),
  variance_test = character(),
  variance_p = numeric(),
  stringsAsFactors = FALSE
)

# ===========================================================
# POSTHOC RESULTS
# ===========================================================
# Pairwise / trend post-hoc comparisons for any significant
# interaction terms, across all variables.
posthoc_table <- data.frame()

# ===========================================================
# POWER STORAGE
# ===========================================================
# Power-analysis results (one row per variable x fixed-effect term,
# and for "perm" variables, per tested effect-size multiplier too).
power_table <- data.frame(
  variable = character(),
  method   = character(),
  effect   = character(),
  effect_multiplier = numeric(),
  power    = numeric(),
  lowerCI  = numeric(),
  upperCI  = numeric(),
  estimate = numeric(),
  nsim     = numeric(),
  stringsAsFactors = FALSE
)

# -----------------------------
# LOAD DATA
# -----------------------------
# One row per trial/timepoint, with repeated measurements per fish
# (hence the random effect on fishID in the models below).
df <- read_csv("combinedWalkDataSend.csv")

# -----------------------------
# FORMAT VARIABLES
# -----------------------------
# Coerce grouping/predictor columns to the types the models expect:
#   - fishID: factor, so it's usable as a mixed-model grouping variable
#   - environment / training: recoded from 0/1 codes to readable factor
#     labels (0 = Aquatic/Untrained, 1 = Terrestrial/Trained)
#   - flowSpeed: numeric (continuous predictor, e.g. treadmill flow speed)
df <- df %>%
  mutate(
    fishID      = factor(fishID),
    environment = factor(environment, levels = c(0,1),
                         labels = c("Aquatic","Terrestrial")),
    training    = factor(training, levels = c(0,1),
                         labels = c("Untrained","Trained")),
    flowSpeed   = as.numeric(flowSpeed)
  )

# ============================================================
# MODEL FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# run_lmm(): fit a linear mixed model for one response variable,
# print its summary/ANOVA, run post-hoc comparisons on any
# significant two-way interaction, record normality/variance
# diagnostics, and append tidy coefficient rows to the global
# results_table. Returns the fitted model object.
#
# Model formula: response ~ (environment + training + flowSpeed)^2
# + (1 | fishID)
#   - "^2" expands to all main effects plus all pairwise
#     interactions among environment, training, and flowSpeed
#     (i.e. environment + training + flowSpeed + environment:training
#     + environment:flowSpeed + training:flowSpeed), but NOT the
#     three-way interaction.
#   - "(1 | fishID)" is a random intercept per fish, accounting for
#     repeated measurements within the same individual.
# ------------------------------------------------------------
run_lmm <- function(df, response_var) {

  form <- as.formula(
    paste(response_var,
          "~ (environment + training + flowSpeed)^2 + (1 | fishID)")
  )

  model <- lmer(form, data = df)

  cat("\n=============================\n")
  cat("LMM:", response_var, "\n")
  cat("=============================\n")

  print(summary(model))
  anova_table <- anova(model)   # Type III (via lmerTest) F-tests with Satterthwaite df
  print(anova_table)

  # -----------------------------
  # POST HOC TESTING
  # -----------------------------

  cat("\nPost hoc testing:\n")

  # Small helper: safely look up a term's p-value from the ANOVA
  # table, regardless of whether lmerTest reports it under
  # "Pr(>F)" or "Pr(>Chisq)". Returns NA if the term isn't present.
  get_p <- function(term) {
    if (term %in% rownames(anova_table)) {

      # Find correct p-value column (Pr(>F) or Pr(>Chisq))
      p_col <- grep("Pr", colnames(anova_table), value = TRUE)[1]

      return(anova_table[term, p_col])

    } else {
      return(NA)
    }
  }

  p_env_tr <- get_p("environment:training")
  p_env_sp <- get_p("environment:flowSpeed")
  p_tr_sp  <- get_p("training:flowSpeed")

  # For each of the three possible two-way interactions, only run
  # (and store) post-hoc comparisons if that interaction was
  # significant (p < .05). Categorical x categorical interactions
  # use emmeans() + pairs() (pairwise comparisons of estimated
  # marginal means); the categorical x continuous (flowSpeed)
  # interactions use emtrends() + pairs() (comparing the slope of
  # flowSpeed across levels of the categorical factor) instead,
  # since a straight pairwise mean comparison doesn't make sense
  # against a continuous predictor.

  # ---- ENVIRONMENT × TRAINING ----
  if (!is.na(p_env_tr) && p_env_tr < 0.05) {

    cat("\n>>> Significant environment × training interaction\n")

    emm <- emmeans(model, ~ environment * training)
    ph <- as.data.frame(pairs(emm, adjust = "tukey"))
    print(ph)

    ph$variable <- response_var
    ph$method <- "LMM"
    ph$interaction <- "environment:training"

    posthoc_table <<- rbind(posthoc_table, ph)
  }

  # ---- ENVIRONMENT × FLOW SPEED ----
  if (!is.na(p_env_sp) && p_env_sp < 0.05) {

    cat("\n>>> Significant environment × flowSpeed interaction\n")

    emm <- emtrends(model, ~ environment, var = "flowSpeed")
    ph <- as.data.frame(pairs(emm))
    print(ph)

    ph$variable <- response_var
    ph$method <- "LMM"
    ph$interaction <- "environment:flowSpeed"

    posthoc_table <<- rbind(posthoc_table, ph)
  }

  # ---- TRAINING × FLOW SPEED ----
  if (!is.na(p_tr_sp) && p_tr_sp < 0.05) {

    cat("\n>>> Significant training × flowSpeed interaction\n")

    emm <- emtrends(model, ~ training, var = "flowSpeed")
    ph <- as.data.frame(pairs(emm))
    print(ph)

    ph$variable <- response_var
    ph$method <- "LMM"
    ph$interaction <- "training:flowSpeed"

    posthoc_table <<- rbind(posthoc_table, ph)
  }

  if (all(is.na(c(p_env_tr, p_env_sp, p_tr_sp)) |
          c(p_env_tr, p_env_sp, p_tr_sp) >= 0.05)) {
    cat("\n>>> No significant interactions — no post hoc needed\n")
  }
  print(r2(model))   # marginal/conditional R^2 for the mixed model (performance::r2)

  # -----------------------------
  # Diagnostics: residual normality (Shapiro-Wilk, via
  # performance::check_normality) and heteroscedasticity
  # (Breusch-Pagan-style test, via check_heteroscedasticity).
  # Both return an object whose printed/formatted text embeds the
  # p-value as "... p = 0.0123 ..."; the regex below pulls that
  # numeric value back out since these functions don't expose a
  # tidy $p.value field directly.
  # -----------------------------
  cat("\nDiagnostics:\n")
  norm <- check_normality(model)
  homo <- check_heteroscedasticity(model)

  print(norm)
  print(homo)

  # Extract p-values safely
  norm_p <- as.numeric(gsub(".*p = ([0-9\\.eE-]+).*", "\\1", format(norm)))
  homo_p <- as.numeric(gsub(".*p = ([0-9\\.eE-]+).*", "\\1", format(homo)))

  diag_row <- data.frame(
    variable = response_var,
    method   = "LMM",
    normality_test = "Shapiro-Wilk",
    normality_p = norm_p,
    variance_test = "Breusch-Pagan",
    variance_p = homo_p
  )

  diagnostics_table <<- rbind(diagnostics_table, diag_row)

  # -----------------------------
  # SAVE LMM RESULTS
  # -----------------------------
  # Pull the fixed-effect coefficient table (estimate, SE, df,
  # t/statistic, p-value) out of the model summary.
  coefs <- summary(model)$coefficients
  coefs_df <- as.data.frame(coefs)

  coefs_df$effect <- rownames(coefs_df)
  # -----------------------------
  # STANDARDIZE EFFECT NAMES (MATCH PERM ANOVA)
  # -----------------------------
  # lmer/lmerTest names dummy-coded terms after their factor level
  # (e.g. "environmentTerrestrial", "trainingTrained"), while the
  # permutation ANOVA (aovperm) reports plain factor names (e.g.
  # "environment", "training"). These gsub() calls rewrite the LMM's
  # coefficient names to match that shared naming convention, so
  # results from both methods can later be compared/filtered by the
  # same `effect` labels.

  coefs_df$effect <- gsub("environmentTerrestrial", "environment", coefs_df$effect)
  coefs_df$effect <- gsub("trainingTrained", "training", coefs_df$effect)

  # These three gsub() calls are no-ops (replacing a pattern with
  # itself) — left in intentionally as a defensive/explicit "these
  # interaction names are already in the target format" marker.
  coefs_df$effect <- gsub("environment:training", "environment:training", coefs_df$effect) # safe
  coefs_df$effect <- gsub("environment:flowSpeed", "environment:flowSpeed", coefs_df$effect)
  coefs_df$effect <- gsub("training:flowSpeed", "training:flowSpeed", coefs_df$effect)

  # Fix interaction naming from LMM style
  # (After the two single-term gsub() calls above already renamed
  # "environmentTerrestrial" -> "environment" etc. within interaction
  # terms too, so by this point terms like
  # "environmentTerrestrial:trainingTrained" have actually already
  # become "environment:training" — these three lines are a second,
  # redundant safety pass in case any raw LMM-style names slipped
  # through unchanged.)
  coefs_df$effect <- gsub("environmentTerrestrial:trainingTrained", "environment:training", coefs_df$effect)
  coefs_df$effect <- gsub("environmentTerrestrial:flowSpeed", "environment:flowSpeed", coefs_df$effect)
  coefs_df$effect <- gsub("trainingTrained:flowSpeed", "training:flowSpeed", coefs_df$effect)
  coefs_df$variable <- response_var
  coefs_df$method <- "LMM"

  colnames(coefs_df) <- c("estimate", "std_error", "df",
                          "statistic", "p_value",
                          "effect", "variable", "method")

  # Append to global table
  results_table <<- rbind(results_table, coefs_df)

  return(model)
}
# ------------------------------------------------------------
# run_log_lmm(): wraps run_lmm() to fit the model on log(y) instead
# of y. Filters out non-positive values first (log is undefined for
# y <= 0), builds a "log_<response_var>" column, fits via run_lmm(),
# then relabels the stored results/diagnostics rows so they appear
# under the *original* variable name with method "logLMM" (rather
# than the literal "log_<var>" name with method "LMM").
# ------------------------------------------------------------
run_log_lmm <- function(df, response_var) {

  log_var <- paste0("log_", response_var)

  df <- df %>%
    filter(.data[[response_var]] > 0)

  df[[log_var]] <- log(df[[response_var]])

  # Run model using log variable
  model <- run_lmm(df, log_var)

  # -----------------------------
  # FIX VARIABLE + METHOD LABELS
  # -----------------------------
  # run_lmm() just wrote rows with variable == log_var and
  # method == "LMM"; rewrite those specific rows in place so the
  # output tables read as variable == response_var, method == "logLMM".
  results_table <<- transform(
    results_table,
    variable = ifelse(variable == log_var, response_var, variable),
    method   = ifelse(variable == response_var & method == "LMM", "logLMM", method)
  )

  diagnostics_table <<- transform(
    diagnostics_table,
    variable = ifelse(variable == log_var, response_var, variable),
    method   = ifelse(variable == response_var & method == "LMM", "logLMM", method)
  )

  return(model)
}
# ------------------------------------------------------------
# run_perm(): fit a permutation-based ANOVA (no random effect term —
# aovperm doesn't support one directly; instead, permutations are
# constrained/clustered by fishID via `cluster =`, which keeps a
# given fish's rows together when shuffling, approximating the
# repeated-measures structure). Used for variables where the LMM's
# parametric assumptions are a poor fit.
# ------------------------------------------------------------
run_perm <- function(df, response_var, n_permutations = 5000) {

  form <- as.formula(
    paste(response_var,
          "~ (environment + training + flowSpeed)^2")
  )

  set.seed(123)

  model <- aovperm(
    form,
    data    = df,
    cluster = df$fishID,
    np      = n_permutations
  )

  cat("\n=============================\n")
  cat("PERMUTATION ANOVA:", response_var, "\n")
  cat("=============================\n")

  perm_summary <- summary(model)
  print(perm_summary)

  # -----------------------------
  # SAVE PERMUTATION RESULTS
  # -----------------------------

  perm_df <- as.data.frame(perm_summary)

  # Remove residual row if present
  perm_df <- perm_df[!grepl("Residuals", rownames(perm_df)), ]

  # Safely grab statistic column
  stat_col <- grep("F", colnames(perm_df), value = TRUE)[1]

  # Prefer resampled p-value
  # (aovperm can report either a "resampled" permutation p-value or
  # a parametric "Pr(>F)"-style column depending on settings; prefer
  # the permutation-based one when available.)
  p_col <- grep("resampled", colnames(perm_df), value = TRUE)
  if (length(p_col) == 0) {
    p_col <- grep("Pr", colnames(perm_df), value = TRUE)
  }
  p_col <- p_col[1]

  perm_df <- perm_df[, c(stat_col, p_col)]
  colnames(perm_df) <- c("statistic", "p_value")

  perm_df$effect <- rownames(perm_df)
  perm_df$variable <- response_var
  perm_df$method <- "perm"

  # -----------------------------
  # STANDARDIZE COLUMN STRUCTURE
  # -----------------------------
  # aovperm results have no coefficient estimate/SE/df in the same
  # sense as an LMM, so these columns are added as NA purely so
  # perm_df has the same column layout as the LMM coefs_df and can
  # be rbind()'d into the shared results_table.
  perm_df$estimate  <- NA
  perm_df$std_error <- NA
  perm_df$df        <- NA

  perm_df <- perm_df[, c(
    "estimate",
    "std_error",
    "df",
    "statistic",
    "p_value",
    "effect",
    "variable",
    "method"
  )]

  # Append
  results_table <<- rbind(results_table, perm_df)

  # -----------------------------
  # SAVE DIAGNOSTICS (NA for perm)
  # -----------------------------
  # Permutation ANOVA makes no distributional assumptions to check,
  # so normality/variance diagnostics are simply recorded as NA —
  # this row exists mainly to keep every analyzed variable
  # represented in diagnostics_table.
  diag_row <- data.frame(
    variable = response_var,
    method   = "perm",
    normality_test = NA,
    normality_p = NA,
    variance_test = NA,
    variance_p = NA
  )

  diagnostics_table <<- rbind(diagnostics_table, diag_row)

  return(model)
}
# ============================================================
# MASTER ANALYSIS PLAN (MANUAL DECISIONS)
# ============================================================
# Which of the three methods (LMM / logLMM / perm) to use for each
# of the 8 outcome variables. These choices were made manually
# (presumably from prior exploratory diagnostics), not derived
# automatically within this script.
analysis_plan <- data.frame(
  variable = c(
    "strideLength_BL",
    "finEffort",
    "tailEffort",
    "totalEffort",
    "tailMeanAmp_BL",
    "finMeanAmp_BL",
    "tailMeanFreq_Hz",
    "finMeanFreq_Hz"
  ),

  method = c(
    "perm",      # stride length
    "LMM",      # finEffort
    "perm",     # tailEffort
    "perm",     # totalEffort
    "LMM",      # tailMeanAmp_BL
    "perm",     # finMeanAmp_BL
    "logLMM",   # tailMeanFreq_Hz
    "logLMM"   # finMeanFreq_Hz
  )
)
# ============================================================
# RUN ALL ANALYSES
# ============================================================
# -----------------------------
# PERMUTATION POWER FUNCTIONS
# -----------------------------

# ------------------------------------------------------------
# simulate_data(): build a synthetic response column
# (`response_sim`) for power simulation. Starting from the real,
# complete-case values of `response_var`, it adds:
#   1. A specified `effect_size` applied through the predictor
#      combination named by `effect_type` (e.g. "environment" adds
#      effect_size for Terrestrial rows only; "env_train" adds it
#      only where BOTH environment==Terrestrial AND training==Trained,
#      approximating an interaction effect).
#   2. A per-fish random intercept (fish_effects), so the simulated
#      data still has fish-level correlation structure like the
#      real data.
#   3. Residual noise scaled off the real data's SD.
# This produces a dataset where the *true* effect size for
# `effect_type` is known exactly (it's injected), so re-running the
# permutation test on many such simulated datasets and checking how
# often it detects p < .05 estimates that test's statistical power.
# ------------------------------------------------------------
simulate_data <- function(df, response_var, effect_size, effect_type) {

  df_sim <- df

  y_base <- df[[response_var]]

  # REMOVE NA rows (critical)
  keep <- !is.na(y_base) &
    !is.na(df$environment) &
    !is.na(df$training) &
    !is.na(df$flowSpeed) &
    !is.na(df$fishID)

  df_sim <- df_sim[keep, ]
  y_base <- y_base[keep]

  sd_val <- sd(y_base)
  if (is.na(sd_val) || sd_val == 0) sd_val <- 1   # guard against a degenerate/constant variable

  # predictors
  env_numeric   <- as.numeric(df_sim$environment == "Terrestrial")
  train_numeric <- as.numeric(df_sim$training == "Trained")
  speed_scaled  <- as.numeric(scale(df_sim$flowSpeed))

  # fish random effect
  # Draw one random intercept per unique fish (SD = 0.3 x the
  # response's overall SD), so simulated observations from the same
  # fish are correlated, similar to the real repeated-measures data.
  fish_ids <- unique(df_sim$fishID)
  fish_effects <- rnorm(length(fish_ids), 0, sd_val * 0.3)
  names(fish_effects) <- fish_ids

  # Build the injected effect according to which term is being
  # power-tested. Main effects scale with a single 0/1 (or
  # standardized continuous) predictor; interaction effects
  # ("env_train", "env_speed", "train_speed") scale with the
  # *product* of the two relevant predictors, i.e. the effect is
  # only present when both conditions hold (or scales with both).
  # Falls back to a zero effect for any unrecognized effect_type.
  effect <- switch(effect_type,
                   "environment"  = effect_size * env_numeric,
                   "training"     = effect_size * train_numeric,
                   "speed"        = effect_size * speed_scaled,
                   "env_train"    = effect_size * env_numeric * train_numeric,
                   "env_speed"    = effect_size * env_numeric * speed_scaled,
                   "train_speed"  = effect_size * train_numeric * speed_scaled,
                   rep(0, length(y_base)))

  # Simulated response = real baseline value + injected effect +
  # per-fish random intercept + residual noise (SD = 0.5 x the
  # response's overall SD).
  df_sim$response_sim <-
    y_base +
    effect +
    fish_effects[as.character(df_sim$fishID)] +
    rnorm(length(y_base), 0, sd_val * 0.5)

  return(df_sim)
}
# ------------------------------------------------------------
# run_perm_test(): fit the permutation ANOVA to one simulated
# dataset (from simulate_data()) and return the p-value for a
# single term of interest. Used inside the power-simulation loop
# below, called many times (once per simulated replicate).
# ------------------------------------------------------------
run_perm_test <- function(df_sim, term) {

  model <- aovperm(
    response_sim ~ (environment + training + flowSpeed)^2,
    data = df_sim,
    np = 5000,
    method = "freedman_lane"
  )

  # Extract table directly
  tab <- model$table

  # Identify correct p-value column
  p_col <- grep("resampled", colnames(tab), value = TRUE)
  if (length(p_col) == 0) p_col <- grep("Pr", colnames(tab), value = TRUE)
  p_col <- p_col[1]

  # Safe extraction
  if (term %in% rownames(tab)) {
    return(as.numeric(tab[term, p_col]))
  } else {
    return(NA)
  }
}
# ============================================================
# MAIN LOOP: run the assigned method for every variable in the
# analysis plan, plus a matching power analysis.
# ============================================================
results <- list()
for (i in 1:nrow(analysis_plan)) {

  var <- analysis_plan$variable[i]
  method <- analysis_plan$method[i]

  cat("\n\n====================================\n")
  cat("Running:", var, "| Method:", method, "\n")
  cat("====================================\n")

  # ================================================
  # METHOD: LMM
  # ================================================
  if (method == "LMM") {

    results[[var]] <- run_lmm(df, var)

    # -----------------------------
    # POWER ANALYSIS
    # -----------------------------
    # Refit the same LMM (with REML = FALSE, as simr recommends for
    # power simulation) on complete cases only, then use simr's
    # powerSim() to estimate power for each fixed-effect term by
    # repeatedly resimulating the response from the fitted model and
    # refitting, at nsim = 1000 simulations per term.

    response_var <- var

    df_clean <- df %>%
      drop_na(all_of(c(response_var, "environment", "training", "flowSpeed", "fishID")))

    df_clean$response_tmp <- df_clean[[response_var]]

    model_power <- lmer(
      response_tmp ~ (environment + training + flowSpeed)^2 + (1 | fishID),
      data = df_clean,
      REML = FALSE
    )


    terms_to_test <- c(
      "environmentTerrestrial",
      "trainingTrained",
      "flowSpeed",
      "environmentTerrestrial:trainingTrained",
      "environmentTerrestrial:flowSpeed",
      "trainingTrained:flowSpeed"
    )

    coefs <- names(fixef(model_power))

    for (term in terms_to_test) {

      cat("\nPower Analysis for:", term, "\n")

      if (!(term %in% coefs)) {
        cat("Skipping", term, "- not in model\n")
        next
      }

      # powerSim() with a z-test on the fixed effect; wrapped in
      # tryCatch so one term's simulation failure doesn't abort the
      # whole loop — a failed simulation just records NA power below.
      pwr <- tryCatch({
        powerSim(
          model_power,
          test = fixed(term, "z"),
          nsim = 1000
        )
      }, error = function(e) {
        cat("PowerSim failed for", term, "\n")
        return(NULL)
      })

      print(pwr)

      # SAFE extraction
      # simr's powerSim() result doesn't expose power/CI as plain
      # numeric fields in a stable way here, so its printed summary
      # (which looks like "Power for predictor... : 45.20% (42.10,
      # 48.35) ...") is captured as text and the first three numbers
      # on the line containing "%" are parsed out as power/lowerCI/
      # upperCI (as percentages, then converted to proportions).
      pwr_out <- if (is.null(pwr)) {
        list(power = NA, lower = NA, upper = NA)
      } else {

        txt <- capture.output(print(pwr))

        power_val <- NA
        lower <- NA
        upper <- NA

        for (line in txt) {
          if (grepl("%", line)) {

            nums <- regmatches(line, gregexpr("[0-9]+\\.?[0-9]*", line))[[1]]

            if (length(nums) >= 3) {
              power_val <- as.numeric(nums[1]) / 100
              lower     <- as.numeric(nums[2]) / 100
              upper     <- as.numeric(nums[3]) / 100
              break
            }
          }
        }

        list(power = power_val, lower = lower, upper = upper)
      }

      # Relabel dummy-coded term names to the shared "environment" /
      # "training" / "environment:training" style used elsewhere.
      effect_clean <- term

      effect_clean <- gsub("environmentTerrestrial", "environment", effect_clean)
      effect_clean <- gsub("trainingTrained", "training", effect_clean)
      effect_clean <- gsub("environmentTerrestrial:trainingTrained", "environment:training", effect_clean)
      effect_clean <- gsub("environmentTerrestrial:flowSpeed", "environment:flowSpeed", effect_clean)
      effect_clean <- gsub("trainingTrained:flowSpeed", "training:flowSpeed", effect_clean)

      power_table <<- rbind(power_table, data.frame(
        variable = response_var,
        method   = "power_LMM",
        effect   = effect_clean,
        effect_multiplier = NA,
        power    = round(pwr_out$power, 3),
        lowerCI  = round(pwr_out$lower, 3),
        upperCI  = round(pwr_out$upper, 3),
        estimate = fixef(model_power)[term],
        nsim     = 1000,
        stringsAsFactors = FALSE
      ))
    }
  }

  # ================================================
  # METHOD: logLMM
  # ================================================
  # Same power-analysis approach as the LMM branch above, but fit on
  # the log-transformed response (mirrors what run_log_lmm() did for
  # the main model fit).
  if (method == "logLMM") {

    results[[var]] <- run_log_lmm(df, var)

    # -----------------------------
    # POWER ANALYSIS (LOG LMM)
    # -----------------------------

    response_var <- paste0("log_", var)

    # 1. Create the log variable
    df_temp <- df %>% filter(.data[[var]] > 0)

    # 2. Assign it to a STATIC name that the model will always use
    df_temp$target_var <- log(df_temp[[var]])

    # 3. Model the static name
    model_power <- lmer(
      target_var ~ (environment + training + flowSpeed)^2 + (1 | fishID),
      data = df_temp,
      REML = FALSE
    )

    terms_to_test <- c(
      "environmentTerrestrial",
      "trainingTrained",
      "flowSpeed",
      "environmentTerrestrial:trainingTrained",
      "environmentTerrestrial:flowSpeed",
      "trainingTrained:flowSpeed"
    )

    coefs <- names(fixef(model_power))

    for (term in terms_to_test) {

      cat("\nPower Analysis for:", term, "\n")

      if (!(term %in% coefs)) {
        cat("Skipping", term, "- not in model\n")
        next
      }

      pwr <- tryCatch({
        powerSim(
          model_power,
          test = fixed(term, "z"),
          nsim = 1000
        )
      }, error = function(e) {
        cat("PowerSim failed for", term, "\n")
        return(NULL)
      })

      print(pwr)

      # SAFE extraction
      pwr_out <- if (is.null(pwr)) {
        list(power = NA, lower = NA, upper = NA)
      } else {

        txt <- capture.output(print(pwr))

        power_val <- NA
        lower <- NA
        upper <- NA

        for (line in txt) {
          if (grepl("%", line)) {

            nums <- regmatches(line, gregexpr("[0-9]+\\.?[0-9]*", line))[[1]]

            if (length(nums) >= 3) {
              power_val <- as.numeric(nums[1]) / 100
              lower     <- as.numeric(nums[2]) / 100
              upper     <- as.numeric(nums[3]) / 100
              break
            }
          }
        }

        list(power = power_val, lower = lower, upper = upper)
      }

      effect_clean <- term

      effect_clean <- gsub("environmentTerrestrial", "environment", effect_clean)
      effect_clean <- gsub("trainingTrained", "training", effect_clean)
      effect_clean <- gsub("environmentTerrestrial:trainingTrained", "environment:training", effect_clean)
      effect_clean <- gsub("environmentTerrestrial:flowSpeed", "environment:flowSpeed", effect_clean)
      effect_clean <- gsub("trainingTrained:flowSpeed", "training:flowSpeed", effect_clean)

      power_table <<- rbind(power_table, data.frame(
        variable = var,
        method   = "power_logLMM",
        effect   = effect_clean,
        effect_multiplier = NA,
        power    = round(pwr_out$power, 3),
        lowerCI  = round(pwr_out$lower, 3),
        upperCI  = round(pwr_out$upper, 3),
        estimate = fixef(model_power)[term],
        nsim     = 1000,
        stringsAsFactors = FALSE
      ))
    }
  }

  # ================================================
  # METHOD: perm
  # ================================================
  # Permutation-ANOVA power analysis is a custom simulate-and-count
  # procedure (permuco/aovperm has no built-in power function like
  # simr does for lmer models):
  #   1. Fit a helper LMM (not saved to results_table) purely to get
  #      realistic effect-size estimates (fixef coefficients) for
  #      each term, since the permutation ANOVA itself doesn't
  #      produce comparable effect-size estimates.
  #   2. For each fixed-effect term, test power at 1x, 2x, and 3x the
  #      observed LMM coefficient as the injected effect size.
  #   3. For each of those effect sizes, simulate 1000 synthetic
  #      datasets (via simulate_data(), each with that exact effect
  #      injected), rerun the permutation ANOVA on each
  #      (run_perm_test()), and compute power as the proportion of
  #      simulations where p < .05, with a binomial exact CI.
  if (method == "perm") {

    results[[var]] <- run_perm(df, var)

    cat("\nRunning permutation power for:", var, "\n")

    # -----------------------------
    # FIT LMM TO GET EFFECT SIZES
    # -----------------------------

    lmm_model <- lmer(
      as.formula(paste(var, "~ (environment + training + flowSpeed)^2 + (1 | fishID)")),
      data = df,
      REML = FALSE
    )

    # 1. Extract coefficients from LMM
    coefs <- fixef(lmm_model)

    # 2. Map shorthand to actual numeric values (Beta)
    # Note: Added env_train here to match your effect_terms
    beta_list <- list(
      environment = coefs["environmentTerrestrial"],
      training    = coefs["trainingTrained"],
      speed       = coefs["flowSpeed"],
      env_train   = coefs["environmentTerrestrial:trainingTrained"],
      env_speed   = coefs["environmentTerrestrial:flowSpeed"],
      train_speed = coefs["trainingTrained:flowSpeed"]
    )

    # 3. Map shorthand to ANOVA table row names
    effect_terms <- list(
      environment = "environment",
      training    = "training",
      speed       = "flowSpeed",
      env_train   = "environment:training",
      env_speed   = "environment:flowSpeed",
      train_speed = "training:flowSpeed"
    )

    # -----------------------------
    # RUN SIMULATION POWER
    # -----------------------------

    for (effect_name in names(beta_list)) {
      beta_obs <- beta_list[[effect_name]]
      if (is.na(beta_obs)) next   # skip if this term wasn't estimable in the LMM

      term <- effect_terms[[effect_name]]

      # Testing 1x, 2x, and 3x the observed effect size
      effect_sizes <- c(1, 2, 3) * beta_obs

      for (eff in effect_sizes) {
        cat(sprintf("Simulating effect: %s at size %.3f\n", effect_name, eff))

        set.seed(123)

        # 1000 simulated datasets at this effect size; for each, run
        # the permutation ANOVA and grab the p-value for this term.
        p_vals <- replicate(1000, {

          # Ensure simulate_data actually adds 'eff' to the 'term' in df
          df_sim <- simulate_data(df, var, eff, effect_name)

          # This function MUST return the p-value for the specific 'term'
          p <- tryCatch(run_perm_test(df_sim, term), error = function(e) NA)
          return(p)
        })

        # Calculate Power and Binomial Confidence Intervals
        # Power = proportion of the 1000 simulated p-values that
        # cleared the p < .05 significance threshold.
        successes <- sum(p_vals < 0.05, na.rm = TRUE)
        total_runs <- sum(!is.na(p_vals))

        power_val <- if (total_runs > 0) successes / total_runs else NA

        # Calculate CI using the Binomial Exact method
        ci <- if(total_runs > 0) binom.test(successes, total_runs)$conf.int else c(NA, NA)

        power_table <<- rbind(power_table, data.frame(
          variable = var,
          method   = "power_perm",
          effect   = term,
          effect_multiplier = round(eff / beta_obs, 2),
          power    = round(power_val, 3),
          lowerCI  = round(ci[1], 3),
          upperCI  = round(ci[2], 3),
          estimate = round(eff, 4),
          nsim     = total_runs,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}
# ============================================================
# END
# ============================================================
cat("\n\nALL ANALYSES COMPLETE\n")
# ============================================================
# SAVE RESULTS
# ============================================================
# Write the four accumulated result tables to CSV, alongside the
# source data, for downstream reporting/plotting.
write.csv(results_table,
          "treadmillStats.csv",
          row.names = FALSE)
cat("\nResults saved successfully.\n")

write.csv(diagnostics_table,
          "treadmillDiagnostics.csv",
          row.names = FALSE)
cat("\nDiagnostics saved successfully.\n")

write.csv(posthoc_table,
          "treadmillPosthoc.csv",
          row.names = FALSE)
cat("\nPosthoc results saved successfully.\n")

write.csv(power_table,
          "treadmillPower.csv",
          row.names = FALSE)
cat("\nPower results saved successfully.\n")
