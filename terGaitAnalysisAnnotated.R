# ============================================================
#  TERRESTRIAL GAIT / KINEMATICS ANALYSIS
# ============================================================
#  What this script does, in plain terms:
#    1. Loads a combined dataset of gait/kinematic variables, where
#       each row is (presumably) one trial/individual, labeled by
#       `group`, `rearedCondition`, and `Training`.
#    2. For each kinematic variable in `vars`:
#         a. Removes within-group outliers (1.5*IQR rule).
#         b. Optionally transforms the variable (log or logit) per
#            the `analysis_plan` table, to better satisfy ANOVA
#            assumptions (normality / homogeneity of variance).
#         c. Runs a 2-way ANOVA: rearedCondition * Training.
#         d. Checks assumptions: Levene's test (equal variance) and
#            Shapiro-Wilk on residuals (normality).
#         e. Computes effect sizes (eta^2, Cohen's f) and observed
#            statistical power for each ANOVA term.
#         f. Runs post-hoc pairwise comparisons (via emmeans) only
#            for terms that were significant (p < .05), with the
#            comparison structure depending on whether the
#            interaction term was significant.
#         g. Generates a QQ plot of model residuals for visual
#            normality inspection.
#    3. Cleans up labels (renames "rearedCondition" -> "environment",
#       "Training" -> "training" for readability) and writes four
#       summary CSVs (ANOVA table, post-hoc comparisons, normality/
#       variance diagnostics, and effect size/power) to disk.
# ============================================================

cat("\f")          # Clear the R console (form-feed character; works in RStudio).
rm(list = ls())     # Remove all objects from the current environment for a clean run.

# ============================================================
#  PACKAGES
# ============================================================
# readr    - fast, consistent CSV I/O (read_csv / write_csv)
# dplyr    - data wrangling (select, filter, mutate, group_by, etc.)
# car      - leveneTest() for homogeneity-of-variance testing
# rstatix  - tidy statistical test helpers (loaded for its ANOVA/stats
#            ecosystem conventions; complements broom/emmeans below)
# broom    - tidy() turns base-R model objects (e.g. aov()) into
#            tidy data frames that are easy to bind_rows() together
# emmeans  - estimated marginal means and pairwise post-hoc contrasts
# nlme     - mixed-effects modeling support (loaded for potential use;
#            not directly called in this script's current logic)
# ggplot2  - QQ plot visualization of model residuals
# pwr      - power analysis (pwr.f2.test for ANOVA-style effect sizes)
library(readr)
library(dplyr)
library(car)
library(rstatix)
library(broom)
library(emmeans)
library(nlme)
library(ggplot2)
library(pwr)

# ------------------------------------------------------------
# Load the combined dataset. Each row should contain the raw
# kinematic measurements plus the grouping columns used below:
# `group`, `rearedCondition`, and `Training`.
# ------------------------------------------------------------
df <- read_csv(
  "combinedTerDataSend.csv"
)

# ============================================================
#  VARIABLES
# ============================================================
# The set of kinematic/gait outcome variables to analyze, one at a
# time, in the loop below. Each will get its own ANOVA, post-hoc
# tests, diagnostics, and effect-size/power estimates.
vars <- c(
  "headElevCycleMean_BL",     # mean head elevation per cycle (body lengths)
  "headElevClipMax_BL",       # max head elevation within a clip (body lengths)
  "curvDorsalCycleMean",      # mean dorsal (body) curvature per cycle
  "curvDorsalClipMax",        # max dorsal curvature within a clip
  "mean_stride_length_BL",    # mean stride length (body lengths)
  "mean_stride_speed_BLs",    # mean stride speed (body lengths/second)
  "mean_dutyFactor",          # mean duty factor (proportion of stride in stance phase)
  "mean_stanceDuration_s",    # mean stance-phase duration (seconds)
  "mean_swingDuration_s",     # mean swing-phase duration (seconds)
  "mean_strideDuration_s",    # mean full stride duration (seconds)
  "meanComSpeed_BLs",         # mean center-of-mass speed (body lengths/second)
  "max_speed_BLs",            # max speed observed (body lengths/second)
  "median_speed_BLs"          # median speed observed (body lengths/second)
)

# ============================================================
#  ANALYSIS PLAN (ALL RAW BY DEFAULT)
# ============================================================
# Maps each variable in `vars` (by position) to a transform method:
#   "raw"   - analyze the variable as measured
#   "log"   - analyze log(y + 1e-6) (right-skewed / always-positive data)
#   "logit" - analyze logit-transformed y (data bounded in [0, 1], e.g.
#             a proportion such as duty factor)
# The default is "raw" for every variable; only variables that failed
# assumption checks (normality/variance) in prior exploratory work were
# switched to "log" or "logit" (see inline flags below).
analysis_plan <- data.frame(
  variable = vars,
  method = c(
    "raw",   # headElevCycleMean_BL
    "log",   # headElevClipMax_BL
    "raw",   # curvDorsalCycleMean
    "raw",   # curvDorsalClipMax
    "raw",   # mean_stride_length_BL
    "raw",   # mean_stride_speed_BLs
    "logit", # mean_dutyFactor  <- ONLY CHANGE (proportion variable, bounded 0-1)
    "raw",   # mean_stanceDuration_s
    "log",   # mean_swingDuration_s
    "raw",   # mean_strideDuration_s
    "raw",   # meanComSpeed_BLs
    "raw",   # max_speed_BLs
    "raw"    # median_speed_BLs
  )
)

# ============================================================
#  OUTLIER REMOVAL
# ============================================================
# Removes outliers within each level of `group_col`, using the
# standard Tukey 1.5*IQR fence: any value below Q1 - 1.5*IQR or
# above Q3 + 1.5*IQR is dropped. Uses tidy-eval ({{ }}) so the
# value and grouping columns can be passed in unquoted.
remove_outliers_by_group <- function(data, value_col, group_col) {
  data %>%
    group_by({{ group_col }}) %>%
    mutate(
      Q1 = quantile({{ value_col }}, 0.25, na.rm = TRUE),
      Q3 = quantile({{ value_col }}, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      keep = {{ value_col }} >= Q1 - 1.5 * IQR &
        {{ value_col }} <= Q3 + 1.5 * IQR
    ) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-Q1, -Q3, -IQR, -keep)   # drop helper columns, keep data tidy
}

# ============================================================
#  STORAGE
# ============================================================
# Empty data frames that accumulate results across all variables
# in the loop below, via bind_rows().
anova_results   <- data.frame()   # tidy ANOVA tables (one block per variable)
posthoc_results <- data.frame()   # emmeans pairwise contrasts (only for sig. effects)
diagnostics     <- data.frame()   # normality (Shapiro) + variance (Levene) checks
effect_size_results <- data.frame()  # eta^2, Cohen's f, observed power per effect

# ============================================================
#  LOOP
# ============================================================
# Iterate over every outcome variable in `vars` and run the full
# analysis pipeline (clean -> transform -> test -> summarize) for
# each one independently.
for (v in vars) {

  cat("\n=============================================\n")
  cat("Variable:", v, "\n")
  cat("=============================================\n")

  # ----------------------------
  # Clean
  # ----------------------------
  # Pull just this variable plus the grouping columns, rename the
  # outcome to a generic `y` (simplifies all downstream formulas),
  # and drop rows with missing values for this variable.
  d <- df %>%
    select(all_of(v), group, rearedCondition, Training) %>%
    rename(y = all_of(v)) %>%
    filter(!is.na(y))

  # Apply the 1.5*IQR outlier filter within each `group`.
  d_clean <- remove_outliers_by_group(d, y, group)

  # ----------------------------
  # RAW / LOG SWITCH
  # ----------------------------
  # Look up this variable's transform method from `analysis_plan`
  # and create `y_used` — the (possibly transformed) column that
  # all subsequent tests/models will actually use.
  method <- analysis_plan$method[analysis_plan$variable == v]

  if (method == "log") {

    # Small epsilon avoids log(0) = -Inf for zero-valued observations.
    d_clean$y_used <- log(d_clean$y + 1e-6)
    method_label <- "ANOVA_log"

  } else if (method == "logit") {

    # Clip to (eps, 1-eps) so the logit transform is defined even if
    # any values are exactly 0 or 1 (e.g. duty factor extremes).
    eps <- 1e-6
    y_clipped <- pmin(pmax(d_clean$y, eps), 1 - eps)
    d_clean$y_used <- log(y_clipped / (1 - y_clipped))
    method_label <- "ANOVA_logit"

  } else {

    # "raw": no transform, analyze the outlier-cleaned value directly.
    d_clean$y_used <- d_clean$y
    method_label <- "ANOVA_raw"
  }

  # ----------------------------
  # Levene
  # ----------------------------
  # Test for homogeneity of variance across the rearedCondition *
  # Training cells. A significant result (p < .05) suggests the
  # equal-variance assumption of ANOVA may be violated.
  lev <- leveneTest(y_used ~ rearedCondition * Training, data = d_clean)

  # ----------------------------
  # ANOVA
  # ----------------------------
  # Two-way factorial ANOVA: main effects of rearedCondition and
  # Training, plus their interaction.
  an <- aov(y_used ~ rearedCondition * Training, data = d_clean)

  # ----------------------------
  # Shapiro
  # ----------------------------
  # Shapiro-Wilk normality test on the ANOVA model residuals. A
  # significant result (p < .05) suggests residuals deviate from
  # normality, which would favor a transform (log/logit) or a
  # non-parametric alternative.
  norm_check <- shapiro.test(residuals(an))

  # ----------------------------
  # Diagnostics
  # ----------------------------
  # Record this variable's normality/variance check results so all
  # variables can be compared side-by-side afterward.
  diagnostics <- bind_rows(
    diagnostics,
    data.frame(
      variable       = v,
      method         = method_label,
      normality_test = "Shapiro-Wilk",
      normality_p    = round(norm_check$p.value, 4),
      variance_test  = "Levene",
      variance_p     = round(lev$`Pr(>F)`[1], 4)
    )
  )

  # ----------------------------
  # Power
  # ----------------------------
  # For each ANOVA term, compute:
  #   eta^2    = SS_effect / SS_total            (proportion of variance explained)
  #   Cohen's f = sqrt(eta^2 / (1 - eta^2))       (effect size for F-tests)
  #   power     = pwr.f2.test(...)$power          (observed/post-hoc power at alpha = .05)
  anova_tab <- anova(an)
  SS_total <- sum(anova_tab$"Sum Sq")

  # Helper: compute eta^2 / Cohen's f / power for one named ANOVA
  # term (e.g. "rearedCondition"), returning NULL if that term isn't
  # present in the model's ANOVA table (defensive — keeps the loop
  # from erroring on models missing a term).
  compute_effects <- function(effect_name, label) {

    if (!(effect_name %in% rownames(anova_tab))) return(NULL)

    SS_effect <- anova_tab[effect_name, "Sum Sq"]
    df_effect <- anova_tab[effect_name, "Df"]
    df_resid  <- anova_tab["Residuals", "Df"]

    eta2 <- SS_effect / SS_total
    f <- sqrt(eta2 / (1 - eta2))

    power <- pwr.f2.test(
      u = df_effect,
      v = df_resid,
      f2 = f^2,
      sig.level = 0.05
    )$power

    data.frame(
      effect = label,
      eta2 = eta2,
      cohen_f = f,
      power_observed = power
    )
  }

  # Compute effect sizes/power for the two main effects and the
  # interaction, using friendlier output labels than the raw
  # ANOVA term names.
  effect_results <- bind_rows(
    compute_effects("rearedCondition", "environment"),
    compute_effects("Training", "training"),
    compute_effects("rearedCondition:Training", "environment:training")
  )

  effect_results$variable <- v
  effect_size_results <- bind_rows(effect_size_results, effect_results)

  # ----------------------------
  # Save ANOVA
  # ----------------------------
  # broom::tidy() converts the aov() object into a data frame
  # (term, df, sumsq, statistic, p.value, ...) for easy stacking
  # across variables.
  an_tab <- broom::tidy(an)
  an_tab$variable <- v
  anova_results <- bind_rows(anova_results, an_tab)

  # ----------------------------
  # POST HOC LOGIC
  # ----------------------------
  # Only run pairwise post-hoc comparisons where there's a
  # significant effect to follow up on (p < .05), and choose the
  # comparison structure based on whether the interaction is
  # significant:
  #   - If rearedCondition:Training IS significant, compare
  #     rearedCondition levels within each Training level (and vice
  #     versa), since the interaction means the simple effects need
  #     to be examined separately rather than pooling across the
  #     other factor.
  #   - If the interaction is NOT significant, fall back to testing
  #     each main effect (rearedCondition, Training) on its own,
  #     but only if that main effect itself was significant.
  if ("rearedCondition:Training" %in% rownames(anova_tab) &&
      anova_tab["rearedCondition:Training", "Pr(>F)"] < 0.05) {

    # Simple effects of rearedCondition within each level of Training.
    emm1 <- emmeans(an, ~ rearedCondition | Training)
    tab1 <- as.data.frame(pairs(emm1))
    tab1$variable <- v
    tab1$effect <- "environment|training"

    # Simple effects of Training within each level of rearedCondition.
    emm2 <- emmeans(an, ~ Training | rearedCondition)
    tab2 <- as.data.frame(pairs(emm2))
    tab2$variable <- v
    tab2$effect <- "training|environment"

    posthoc_results <- bind_rows(posthoc_results, tab1, tab2)

  } else {

    # No significant interaction: test main effects independently.
    if ("rearedCondition" %in% rownames(anova_tab) &&
        anova_tab["rearedCondition", "Pr(>F)"] < 0.05) {

      emm_env <- emmeans(an, ~ rearedCondition)
      env_tab <- as.data.frame(pairs(emm_env))
      env_tab$variable <- v
      env_tab$effect <- "environment"

      posthoc_results <- bind_rows(posthoc_results, env_tab)
    }

    if ("Training" %in% rownames(anova_tab) &&
        anova_tab["Training", "Pr(>F)"] < 0.05) {

      emm_train <- emmeans(an, ~ Training)
      train_tab <- as.data.frame(pairs(emm_train))
      train_tab$variable <- v
      train_tab$effect <- "training"

      posthoc_results <- bind_rows(posthoc_results, train_tab)
    }
  }

  # ----------------------------
  # QQ Plot
  # ----------------------------
  # Visual check of residual normality to complement the
  # Shapiro-Wilk test above. Printed to the active graphics device
  # for each variable in turn (e.g. saved into a PDF if the script
  # is run inside a pdf()/dev.off() block, or shown interactively).
  p_qq <- ggplot(data.frame(res = residuals(an)), aes(sample = res)) +
    stat_qq() +
    stat_qq_line(color = "red") +
    ggtitle(paste("QQ:", v)) +
    theme_classic()

  print(p_qq)
}

# ============================================================
#  CLEAN OUTPUT
# ============================================================
# Rename columns and relabel factor names for readability in the
# exported CSVs: "term"/"statistic" -> "effect"/"F", and
# "rearedCondition"/"Training" -> "environment"/"training" wherever
# they appear inside effect-name strings (e.g.
# "rearedCondition:Training" -> "environment:training").
anova_results <- anova_results %>%
  rename(effect = term, F = statistic, p = p.value)
anova_results$effect <- gsub("rearedCondition", "environment", anova_results$effect)
anova_results$effect <- gsub("Training", "training", anova_results$effect)

posthoc_results <- posthoc_results %>%
  rename(p = p.value)

effect_size_results <- effect_size_results %>%
  mutate(
    eta2 = round(eta2, 4),
    cohen_f = round(cohen_f, 4),
    power_observed = round(power_observed, 4)
  )

print(anova_results)

# ============================================================
#  SAVE
# ============================================================
# Write the four accumulated result tables out as CSVs, alongside
# the source data, for downstream reporting/plotting.
out_dir <- "dataOutputPolyp"

write_csv(anova_results,   file.path(out_dir, "terStats.csv"))       # ANOVA main/interaction effects for all variables
write_csv(posthoc_results, file.path(out_dir, "terPosthoc.csv"))     # pairwise post-hoc comparisons (sig. effects only)
write_csv(diagnostics,     file.path(out_dir, "terDiagnostics.csv")) # Shapiro-Wilk + Levene results per variable
write_csv(effect_size_results, file.path(out_dir, "terPower.csv"))   # eta^2, Cohen's f, observed power per effect
