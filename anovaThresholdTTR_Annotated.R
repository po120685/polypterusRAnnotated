# =========================================================
# TIME TO RECOVERY (TTR) ANALYSIS — THRESHOLD METHOD
# TWO-WAY ANOVA VERSION (Treatment x Exercise)
# =========================================================
# PURPOSE
#   Estimate how long each fish takes to return to near-resting
#   oxygen consumption (MO2) after exercise, then test whether
#   recovery time depends on rearing environment (Treatment),
#   exercise training (Exercise), or their interaction.
#
# RECOVERY DEFINITION
#   A fish is "recovered" at the first minute where:
#
#       MO2 <= resting_MO2 * threshold
#
#   AND every later minute in the window also stays at or
#   below that threshold (sustained recovery, not a single dip).
#
# STATISTICAL DESIGN
#   2 x 2 factorial:
#     Treatment: terrestrial vs aquatic (rearing environment)
#     Exercise : ramp vs noRamp         (trained vs untrained)
#   Model: response ~ Treatment * Exercise
#   Tests: Treatment main effect, Exercise main effect,
#          Treatment x Exercise interaction
#   Sums of squares: Type III (robust to unequal group sizes)
#
# OUTPUTS
#   - One TTR value per fish
#   - Sample size per cell
#   - Two-way ANOVA table (Type III)
#   - Assumption checks (normality, equal variance)
#   - Estimated marginal means + Tukey-adjusted comparisons
#   - Boxplot and interaction plot
#
# CAVEATS
#   - Fish that never recover within the window get TTR = NA
#     and are DROPPED before the ANOVA. These are the slowest
#     recoverers, so excluding them can bias group means.
#   - A fish whose only below-threshold reading is the final
#     minute still counts as recovered at that minute.
# =========================================================

cat("\f")       # Clear the R console (RStudio)
rm(list = ls()) # Remove all objects for a clean run

# =========================================================
# LOAD PACKAGES
# =========================================================

library(tidyverse) # Data wrangling (dplyr) and plotting (ggplot2)
library(car)       # Anova() for Type III SS; leveneTest()
library(emmeans)   # Estimated marginal means + post hoc tests

# =========================================================
# USER SETTINGS
# =========================================================

# Input file: one row per fish per minute of recovery.
# Required columns: Tank_ID, Mode, Treatment, Exercise,
# minute, mo2, resting_mo2.
data_file <- "recoverySpreadsheet1min.csv"

# ---------------------------------
# MODE TO ANALYZE
# ---------------------------------
# Run once per locomotor mode by toggling these lines.

#selected_mode <- "walking"
selected_mode <- "swimming"

# ---------------------------------
# RECOVERY THRESHOLD
# ---------------------------------
# Multiplier applied to resting MO2.
#   1.05 = within 5% of resting MO2
#   1.10 = within 10% of resting MO2  (current setting)

threshold <- 1.10

# ---------------------------------
# RECOVERY WINDOW (minutes)
# ---------------------------------
# TTR cannot exceed this value; fish not recovered by
# this minute get TTR = NA.

max_minute <- 30

# ---------------------------------
# LOG TRANSFORM?
# ---------------------------------
# TRUE  = analyze log(TTR)
# FALSE = analyze raw TTR (minutes)
# Default reproduces the earlier analysis: log for swimming,
# raw for walking. For a consistent analysis across modes,
# set this to a single TRUE or FALSE.

use_log <- selected_mode == "swimming"

# =========================================================
# GROUP COLORS
# =========================================================
# Terrestrial-reared = browns; aquatic-reared = blues.

grpColors <- c(
  "TT" = "#8B4513",   # dark brown  (terrestrial, ramp)
  "TU" = "#CC843F",   # tan         (terrestrial, noRamp)
  "AT" = "#1A4C99",   # dark blue   (aquatic, ramp)
  "AU" = "#4DB3FF"    # sky blue    (aquatic, noRamp)
)

# =========================================================
# CONTRASTS FOR TYPE III TESTS
# =========================================================
# Type III sums of squares are only meaningful with
# sum-to-zero ("effects") contrasts. R's default
# (treatment contrasts) gives misleading Type III results.

options(contrasts = c("contr.sum", "contr.poly"))

# =========================================================
# LOAD AND FILTER DATA
# =========================================================

df <- read.csv(data_file)

df <- df %>%
  filter(Mode == selected_mode,   # keep selected locomotor mode
         minute <= max_minute)    # restrict to recovery window

# =========================================================
# CLEAN DATA
# =========================================================
# Convert design variables to factors with a fixed level order
# so tables and plots are consistent.

df$Tank_ID   <- factor(df$Tank_ID)
df$Treatment <- factor(df$Treatment, levels = c("terrestrial", "aquatic"))
df$Exercise  <- factor(df$Exercise,  levels = c("ramp", "noRamp"))

# =========================================================
# CREATE GROUP VARIABLE (for plotting / labels)
# =========================================================
# First letter  = rearing (T = terrestrial, A = aquatic)
# Second letter = exercise (T = trained/ramp, U = untrained/noRamp)
# The ANOVA uses Treatment and Exercise directly; 'group' is
# just a convenient label for plots and tables.

df <- df %>%
  mutate(
    group = case_when(
      Treatment == "terrestrial" & Exercise == "ramp"   ~ "TT",
      Treatment == "terrestrial" & Exercise == "noRamp" ~ "TU",
      Treatment == "aquatic"     & Exercise == "ramp"   ~ "AT",
      Treatment == "aquatic"     & Exercise == "noRamp" ~ "AU"
    ),
    group = factor(group, levels = c("TT", "TU", "AT", "AU"))
  )

# =========================================================
# FLAG RECOVERY AT EACH MINUTE
# =========================================================
# Each fish's cutoff is based on its OWN resting MO2.
# recovered = TRUE when MO2 is at/below the cutoff.

df <- df %>%
  mutate(
    recovery_threshold = resting_mo2 * threshold,
    recovered          = mo2 <= recovery_threshold
  )

# =========================================================
# CALCULATE TTR PER FISH
# =========================================================
# For each fish, step through minutes in order and find the
# FIRST minute from which all remaining minutes are recovered.
# If none, TTR = NA (did not recover within the window).

TTR_results <- list()

for(f in unique(df$Tank_ID)) {
  
  # This fish's data, sorted by time
  fish_df <- df %>%
    filter(Tank_ID == f) %>%
    arrange(minute)
  
  recovered_vec <- fish_df$recovered
  TTR <- NA
  
  for(i in seq_along(recovered_vec)) {
    
    # Sustained recovery: this minute AND all later minutes
    # must be below threshold.
    if(all(recovered_vec[i:length(recovered_vec)])) {
      TTR <- fish_df$minute[i]
      break
    }
  }
  
  # Keep identifying info + TTR
  TTR_results[[as.character(f)]] <- fish_df[1, ] %>%
    select(Tank_ID, Treatment, Exercise, group) %>%
    mutate(TTR = TTR)
}

TTR_all <- bind_rows(TTR_results)

print(TTR_all)  # NA = did not recover within the window

# =========================================================
# RECOVERY SUMMARY (BEFORE DROPPING NAs)
# =========================================================
# How many fish per group recovered vs not. Useful context
# because non-recovered fish are excluded from the ANOVA.

recovery_summary <- TTR_all %>%
  group_by(group) %>%
  summarise(
    n_total     = n(),
    n_recovered = sum(!is.na(TTR)),
    n_not_recov = sum(is.na(TTR)),
    .groups = "drop"
  )

cat("\n--- Recovery summary ---\n")
print(recovery_summary)

# =========================================================
# REMOVE NON-RECOVERED FISH
# =========================================================

TTR_df <- TTR_all %>%
  filter(!is.na(TTR))

# Sample size per cell actually entering the ANOVA.
# Every cell needs n >= 2 to estimate the interaction.
cat("\n--- Cell sizes in ANOVA ---\n")
print(table(TTR_df$Treatment, TTR_df$Exercise))

# =========================================================
# RESPONSE VARIABLE
# =========================================================

if(use_log) {
  cat("\nAnalyzing log(TTR)\n")
  TTR_df$response <- log(TTR_df$TTR)
} else {
  cat("\nAnalyzing raw TTR (minutes)\n")
  TTR_df$response <- TTR_df$TTR
}

# =========================================================
# TWO-WAY ANOVA
# =========================================================
# Fit the full factorial model: main effects of Treatment
# and Exercise plus their interaction.

anova_model <- aov(
  response ~ Treatment * Exercise,
  data = TTR_df
)

# Type III sums of squares via car::Anova().
# Each effect is tested after adjusting for all others,
# so results don't depend on the order of terms, which
# matters here because group sizes are unequal.
# (Base summary(anova_model) would give order-dependent
# Type I results, so it is not used.)
cat("\n--- Two-way ANOVA (Type III SS) ---\n")
anova_table <- Anova(anova_model, type = 3)
print(anova_table)

# HOW TO READ IT
#   Treatment          : do terrestrial vs aquatic fish differ
#                        in TTR (averaged over exercise)?
#   Exercise           : do ramp vs noRamp fish differ
#                        (averaged over rearing)?
#   Treatment:Exercise : does the exercise effect depend on
#                        rearing environment?
#   Ignore the (Intercept) row: it only tests whether mean
#   TTR differs from zero.
#   Report each effect as F(df_effect, df_residual), p.

# =========================================================
# ASSUMPTION TESTS
# =========================================================

# ---------------------------------
# NORMALITY OF RESIDUALS
# ---------------------------------
# Shapiro–Wilk: p > 0.05 -> no evidence of non-normality.
# Low power with small n, so also inspect the Q-Q plot.

cat("\n--- Shapiro-Wilk on residuals ---\n")
print(shapiro.test(residuals(anova_model)))

qqnorm(residuals(anova_model), main = "Q-Q plot of residuals")
qqline(residuals(anova_model))

# ---------------------------------
# HOMOGENEITY OF VARIANCE
# ---------------------------------
# Levene's test (median-centered) across the four cells.
# p > 0.05 -> no evidence variances differ.

cat("\n--- Levene's test ---\n")
print(leveneTest(response ~ Treatment * Exercise, data = TTR_df))

# =========================================================
# POST HOC TESTS
# =========================================================
# Estimated marginal means (model-based group means).
# If use_log = TRUE, type = "response" back-transforms the
# means to minutes (geometric means); contrasts become ratios.

emm_type <- if(use_log) "response" else "link"

# ---------------------------------
# MAIN EFFECTS
# ---------------------------------
# Only interpretable on their own if the interaction is
# NOT significant. Each is averaged over the other factor.

cat("\n--- Treatment main effect ---\n")
print(emmeans(anova_model, pairwise ~ Treatment, type = emm_type))

cat("\n--- Exercise main effect ---\n")
print(emmeans(anova_model, pairwise ~ Exercise, type = emm_type))

# ---------------------------------
# SIMPLE EFFECTS
# ---------------------------------
# Use these if the interaction IS significant: the effect
# of exercise within each rearing environment, and vice versa.

cat("\n--- Exercise effect within each Treatment ---\n")
print(emmeans(anova_model, pairwise ~ Exercise | Treatment, type = emm_type))

cat("\n--- Treatment effect within each Exercise ---\n")
print(emmeans(anova_model, pairwise ~ Treatment | Exercise, type = emm_type))

# ---------------------------------
# ALL FOUR CELLS
# ---------------------------------
# All 6 pairwise comparisons among TT, TU, AT, AU with
# Tukey adjustment (equivalent to the earlier one-way post hoc).

cat("\n--- All pairwise cell comparisons (Tukey) ---\n")
print(emmeans(anova_model, pairwise ~ Treatment * Exercise, type = emm_type))

# =========================================================
# SAVE RESULTS (optional)
# =========================================================

#write.csv(TTR_all,
#          paste0("TTR_threshold_", selected_mode, ".csv"),
#          row.names = FALSE)

# =========================================================
# BOXPLOT
# =========================================================
# Raw TTR (minutes) by group with individual fish overlaid.

p_box <- ggplot(TTR_df,
                aes(x = group, y = TTR, fill = group)) +
  
  # Median and IQR per group; outliers hidden because
  # every fish is plotted as a point
  geom_boxplot(width = 0.7, alpha = 0.85, outlier.shape = NA) +
  
  # Individual fish
  geom_jitter(width = 0.12, height = 0, size = 2.5, alpha = 0.8) +
  
  scale_fill_manual(values = grpColors) +
  
  theme_classic(base_size = 14) +
  
  labs(title = paste("Time To Recovery (TTR) |", selected_mode),
       x = "", y = "TTR (minutes)") +
  
  theme(legend.position = "none")

print(p_box)

# =========================================================
# INTERACTION PLOT
# =========================================================
# Mean TTR (± SE) for each Treatment x Exercise cell.
# Parallel lines = no interaction; crossing or diverging
# lines = the exercise effect differs by rearing environment.

interaction_summary <- TTR_df %>%
  group_by(Treatment, Exercise) %>%
  summarise(
    mean_TTR = mean(TTR),
    se_TTR   = sd(TTR) / sqrt(n()),
    .groups  = "drop"
  )

p_int <- ggplot(interaction_summary,
                aes(x = Exercise, y = mean_TTR,
                    color = Treatment, group = Treatment)) +
  
  geom_point(size = 3, position = position_dodge(width = 0.15)) +
  geom_line(position = position_dodge(width = 0.15)) +
  geom_errorbar(aes(ymin = mean_TTR - se_TTR,
                    ymax = mean_TTR + se_TTR),
                width = 0.1,
                position = position_dodge(width = 0.15)) +
  
  scale_color_manual(values = c("terrestrial" = "#8B4513",
                                "aquatic"     = "#1A4C99")) +
  
  theme_classic(base_size = 14) +
  
  labs(title = paste("Treatment x Exercise |", selected_mode),
       x = "Exercise", y = "Mean TTR (minutes ± SE)",
       color = "Rearing")

print(p_int)