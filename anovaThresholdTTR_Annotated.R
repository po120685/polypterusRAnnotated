# =========================================================
# TIME TO RECOVERY (TTR) ANALYSIS — THRESHOLD METHOD
# =========================================================
# PURPOSE
#   Estimate how long each fish takes to return to near-resting
#   oxygen consumption (MO2) after exercise, then compare that
#   recovery time among four treatment groups.
#
# RECOVERY DEFINITION
#   A fish is "recovered" at the first minute where:
#
#       MO2 <= resting_MO2 * threshold
#
#   AND every later minute in the window also stays at or
#   below that threshold (i.e., recovery must be sustained,
#   not a single dip).
#
# OUTPUTS
#   - One TTR value per fish
#   - One-way ANOVA comparing TTR among groups
#   - Assumption checks (normality, equal variance)
#   - Tukey-adjusted pairwise comparisons
#   - Boxplot of TTR by group
#
# NOTES / CAVEATS
#   - Fish that never meet the criterion within the window
#     get TTR = NA and are DROPPED before the ANOVA. Because
#     these are the slowest recoverers, excluding them can
#     bias group means; a survival analysis (treating them
#     as censored at 30 min) avoids this.
#   - A fish whose only below-threshold reading is the final
#     minute (30) still counts as recovered at minute 30,
#     because there are no later points left to check.
# =========================================================

cat("\f")       # Clear the R console (RStudio)
rm(list = ls()) # Remove all objects from the workspace for a clean run

# =========================================================
# LOAD PACKAGES
# =========================================================

library(tidyverse) # Data wrangling (dplyr) and plotting (ggplot2)
library(car)       # leveneTest() for homogeneity of variance
library(emmeans)   # Estimated marginal means + Tukey post hoc tests

# =========================================================
# USER SETTINGS
# =========================================================

# Input file: one row per fish per minute of recovery.
# Expected columns include: Tank_ID, Mode, Treatment, Exercise,
# minute, mo2, resting_mo2.
data_file <- "recoverySpreadsheet1min.csv"

# ---------------------------------
# MODE TO ANALYZE
# ---------------------------------
# Run the script once per locomotor mode by toggling these lines.

selected_mode <- "walking"
#selected_mode <- "swimming"

# ---------------------------------
# RECOVERY THRESHOLD
# ---------------------------------
# Multiplier applied to resting MO2.
#   1.05 = within 5% of resting MO2
#   1.10 = within 10% of resting MO2  (current setting)

threshold <- 1.10

# =========================================================
# GROUP COLORS
# =========================================================
# Fixed colors so groups look the same across all figures.
# Terrestrial-reared groups = browns; aquatic-reared = blues.

grpColors <- c(
  "TT" = "#8B4513",   # dark brown  (terrestrial, ramp)
  "TU" = "#CC843F",   # tan         (terrestrial, noRamp)
  "AT" = "#1A4C99",   # dark blue   (aquatic, ramp)
  "AU" = "#4DB3FF"    # sky blue    (aquatic, noRamp)
)

# =========================================================
# LOAD DATA
# =========================================================

df <- read.csv(data_file)

# =========================================================
# FILTER MODE
# =========================================================
# Keep only rows for the locomotor mode chosen above.

df <- df %>%
  filter(Mode == selected_mode)

# =========================================================
# LIMIT TO FIRST 30 MINUTES
# =========================================================
# Restrict the recovery window to minutes 0–30. TTR can
# therefore never exceed 30, and fish that haven't recovered
# by minute 30 will be NA.

df <- df %>%
  filter(minute <= 30)

# =========================================================
# CLEAN DATA
# =========================================================
# Convert ID and design variables to factors so R treats
# them as categories rather than text/numbers.

df$Tank_ID   <- factor(df$Tank_ID)
df$Treatment <- factor(df$Treatment)
df$Exercise  <- factor(df$Exercise)

# =========================================================
# CREATE GROUP VARIABLE
# =========================================================
# Combine the two design factors (Treatment x Exercise) into
# a single 4-level grouping variable:
#   First letter  = rearing environment (T = terrestrial, A = aquatic)
#   Second letter = exercise (T = trained/ramp, U = untrained/noRamp)

df$group <- case_when(
  
  df$Treatment == "terrestrial" &
    df$Exercise == "ramp" ~ "TT",
  
  df$Treatment == "terrestrial" &
    df$Exercise == "noRamp" ~ "TU",
  
  df$Treatment == "aquatic" &
    df$Exercise == "ramp" ~ "AT",
  
  df$Treatment == "aquatic" &
    df$Exercise == "noRamp" ~ "AU"
)

# Set a fixed level order so tables and plots always list
# groups as TT, TU, AT, AU.
df$group <- factor(
  df$group,
  levels = c("TT", "TU", "AT", "AU")
)

# =========================================================
# CALCULATE THRESHOLD
# =========================================================
# Each fish's recovery cutoff is based on its OWN resting
# MO2, so the criterion is individual-specific.

df <- df %>%
  
  mutate(
    recovery_threshold = resting_mo2 * threshold
  )

# =========================================================
# DETERMINE RECOVERY STATUS
# =========================================================
# For every minute, flag whether MO2 is at/below the cutoff.
# TRUE = below threshold at that minute; FALSE = still elevated.

df <- df %>%
  
  mutate(
    recovered = mo2 <= recovery_threshold
  )

# =========================================================
# CALCULATE TTR PER FISH
# =========================================================
# For each fish, step through minutes in order and find the
# FIRST minute from which all remaining minutes are recovered.
# That minute is the fish's TTR. If no such minute exists,
# TTR stays NA (never recovered within 30 min).

TTR_results <- list()            # Container for one row per fish

fish_ids <- unique(df$Tank_ID)   # Each Tank_ID = one fish

for(f in fish_ids) {
  
  # Pull this fish's data, sorted by time
  fish_df <- df %>%
    
    filter(Tank_ID == f) %>%
    
    arrange(minute)
  
  # Vector of TRUE/FALSE recovery flags across minutes
  recovered_vec <- fish_df$recovered
  
  # Default: not recovered
  TTR <- NA
  
  for(i in 1:length(recovered_vec)) {
    
    # Recovery flags from minute i to the end of the window
    current_and_remaining <-
      recovered_vec[i:length(recovered_vec)]
    
    # Recovery only counts if ALL remaining points stay
    # recovered (sustained recovery). The first i where this
    # holds is the TTR; stop searching once found.
    # NOTE: at the last minute only one point remains, so a
    # single low final reading is enough to give TTR = 30.
    
    if(all(current_and_remaining == TRUE)) {
      
      TTR <- fish_df$minute[i]
      
      break
    }
  }
  
  # Keep the fish's first row (for its ID/group info)
  # and attach the TTR value
  temp <- fish_df[1, ]
  
  temp$TTR <- TTR
  
  TTR_results[[as.character(f)]] <- temp
}

# =========================================================
# COMBINE TTR TABLE
# =========================================================
# Stack the per-fish rows into one data frame.

TTR_df <- bind_rows(TTR_results)

# Keep only identifying columns and the TTR result
TTR_df <- TTR_df %>%
  
  select(
    Tank_ID,
    Treatment,
    Exercise,
    group,
    TTR
  )

print(TTR_df)  # Inspect: NA = did not recover within 30 min

# =========================================================
# REMOVE NAs
# =========================================================
# Drop fish that never recovered. The ANOVA below is
# therefore run ONLY on fish that recovered within the window.
# (See caveat in header: this excludes the slowest fish.)

TTR_df <- TTR_df %>%
  
  filter(!is.na(TTR))

# =========================================================
# CONDITIONAL LOG TRANSFORM
# =========================================================
# Swimming TTR is log-transformed (intended to improve
# normality of residuals); walking TTR is analyzed raw.
# NOTE: if reporting both modes, consider using the same
# scale for both so the analyses are directly comparable.

if(selected_mode == "swimming") {
  
  cat("\nApplying log transform to swimming TTR\n")
  
  TTR_df$response <- log(TTR_df$TTR)
  
} else {
  
  cat("\nUsing raw TTR for walking\n")
  
  TTR_df$response <- TTR_df$TTR
}

# =========================================================
# ONE-WAY ANOVA
# =========================================================
# Tests whether mean TTR differs among the four groups.
# NOTE: this is a ONE-way ANOVA on the combined 'group'
# factor. To test rearing and exercise effects separately
# (and their interaction), use instead:
#     aov(response ~ Treatment * Exercise, data = TTR_df)

anova_model <- aov(
  response ~ group,
  data = TTR_df
)

summary(anova_model)  # F statistic and p-value for the group effect

# =========================================================
# ASSUMPTION TESTS
# =========================================================

# ---------------------------------
# NORMALITY
# ---------------------------------
# Shapiro–Wilk test on the model residuals.
# p > 0.05 -> no evidence residuals depart from normality.
# (Low power with small samples; also check a Q-Q plot.)

shapiro_result <- shapiro.test(
  residuals(anova_model)
)

print(shapiro_result)

# ---------------------------------
# HOMOGENEITY OF VARIANCE
# ---------------------------------
# Levene's test (median-centered) for equal variance
# across groups.
# p > 0.05 -> no evidence variances differ among groups.

levene_result <- leveneTest(
  response ~ group,
  data = TTR_df
)

print(levene_result)

# =========================================================
# POST HOC TESTS
# =========================================================
# Estimated marginal (model-based) means for each group,
# plus all pairwise group differences with Tukey adjustment
# for multiple comparisons.
# For the swimming (log) model, estimates are on the log scale;
# add type = "response" to emmeans() to back-transform.

emmeans_results <- emmeans(
  anova_model,
  pairwise ~ group
)

print(emmeans_results)

# =========================================================
# SAVE RESULTS
# =========================================================
# Uncomment to export the per-fish TTR table, e.g.
# "TTR_threshold_walking.csv".

#write.csv(
#TTR_df,
#paste0(
#"TTR_threshold_",
#selected_mode,
#".csv"
#),
#row.names = FALSE
#)

# =========================================================
# BOXPLOT
# =========================================================
# Raw TTR (minutes, untransformed) by group, with individual
# fish overlaid as jittered points.

ggplot(TTR_df,
       aes(x = group,
           y = TTR,
           fill = group)) +
  
  # Box = median and interquartile range per group.
  # Outliers hidden here because every fish is shown as a point.
  geom_boxplot(
    width = 0.7,
    alpha = 0.85,
    outlier.shape = NA
  ) +
  
  # Individual fish, spread horizontally so points don't overlap
  geom_jitter(
    width = 0.12,
    size = 2.5,
    alpha = 0.8
  ) +
  
  # Apply the custom group colors to box fills
  scale_fill_manual(values = grpColors) +
  
  # NOTE: no 'color' aesthetic is mapped, so this line does
  # nothing and triggers a harmless warning. Delete it, or add
  # aes(color = group) to geom_jitter() to color the points.
  scale_color_manual(values = grpColors) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = paste(
      "Time To Recovery (TTR)",
      "|",
      selected_mode
    ),
    x = "",
    y = "TTR (minutes)"
  ) +
  
  # Legend is redundant since groups are labeled on the x-axis
  theme(
    legend.position = "none"
  )