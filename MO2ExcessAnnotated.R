# =========================================================
# TOTAL ENERGY ANALYSIS
# SEPARATE 2-WAY ANOVAS FOR EACH MODE
# =========================================================
#
# PURPOSE OF THIS SCRIPT
# -----------------------------------------------------------------------
# For each locomotion Mode ("swimming" and "walking") separately, this
# script:
#   1. Runs a classic 2-way ANOVA (Treatment x Exercise) on total energy.
#   2. Runs a permutation-based 2-way ANOVA as a robust alternative that
#      doesn't rely on normality assumptions.
#   3. Checks ANOVA assumptions (normality of residuals, equal variance).
#   4. Runs post-hoc pairwise comparisons (emmeans).
#   5. Makes boxplots of total energy by group, split by Mode.
#   6. Saves the permutation ANOVA results tables to CSV.
#   7. Runs a simulation-based power analysis for the Treatment:Exercise
#      interaction term, across a range of hypothetical effect sizes, so
#      you can see how much data/effect size you'd need to reliably
#      detect that interaction with the permutation ANOVA.
#
# =========================================================

# cat("\f") sends a "form feed" character to the console, which clears
# the RStudio console pane (cosmetic only -- doesn't affect any objects).
cat("\f")

# Remove every object currently in the global environment so the script
# runs from a clean slate and old variables can't accidentally leak in.
rm(list = ls())

# -----------------------------
# LOAD PACKAGES
# -----------------------------
library(tidyverse)   # data wrangling (dplyr, etc.) + ggplot2 for plotting
library(car)         # provides leveneTest() for homogeneity-of-variance testing
library(emmeans)     # "estimated marginal means" -- post-hoc/pairwise comparisons
library(permuco)     # permutation-based ANOVA (aovperm), robust to non-normal data
library(svglite)     # lets ggsave() write vector .svg files

# -----------------------------
# LOAD DATA
# -----------------------------

# Read in the per-fish summed metabolic-rate/energy data 
df <- read.csv("MO2TTRSumPerFish.csv")

# Folder where all output figures/tables from this script will be saved.
out_dir <- "dataOutputPolyp"

# -----------------------------
# FORMAT VARIABLES
# -----------------------------

# Tank_ID is an identifier, not a continuous number, so make it a factor
# (categorical) rather than leaving it as an integer/character.
df$Tank_ID <- factor(df$Tank_ID)

# Mode = the locomotion type being measured (swimming vs. walking trial).
# Setting explicit `levels` fixes the category order (rather than relying
# on default alphabetical order), which matters for how results/plots are
# ordered and how contrasts are coded internally.
df$Mode <- factor(df$Mode,
                  levels = c("swimming", "walking"))

# Treatment = the rearing/acclimation condition (aquatic vs. terrestrial).
df$Treatment <- factor(df$Treatment,
                       levels = c("aquatic", "terrestrial"))

# Exercise = whether the fish underwent the training or not ramp = trained, noRamp = untrained.
df$Exercise <- factor(df$Exercise,
                      levels = c("noRamp", "ramp"))

# =========================================================
# LOG TRANSFORM
# (recommended for energy data)
# =========================================================

# A log transformation of the data
df$log_energy <- log(df$total_energy_J)

# =========================================================
# SPLIT DATASETS
# =========================================================

# Because the question is "does Treatment/Exercise affect energy use
# within swimming, and within walking" 
df_swim <- subset(df, Mode == "swimming")

df_walk <- subset(df, Mode == "walking")

# =========================================================
# SWIMMING 2-WAY ANOVA
# =========================================================

# Standard (OLS) 2-way ANOVA: tests main effects of Treatment and
# Exercise, plus their interaction, on total_energy_J within the
# swimming subset.
model_swim <- aov(total_energy_J ~ Treatment * Exercise,
                  data = df_swim)

# Print the ANOVA table (F-statistics, p-values) for each term.
summary(model_swim)

# =========================================================
# SWIMMING
# =========================================================

# Permutation ANOVA: instead of relying on the F-distribution (which
# assumes normally distributed residuals), this repeatedly shuffles
# (permutes) the data np = 5000 times to build an empirical null
# distribution for each F-statistic. More robust when normality is in
# doubt, at the cost of being computationally heavier.
perm_swim <- aovperm(
  total_energy_J ~ Treatment * Exercise,
  data = df_swim,
  np = 5000
)

summary(perm_swim)


# -----------------------------
# ASSUMPTION TESTS
# -----------------------------

# Normality: Shapiro-Wilk test on the residuals of the parametric ANOVA
# model. A significant (low p-value) result suggests residuals deviate
# from normality, motivating use of the permutation ANOVA above.
shapiro.test(residuals(model_swim))

# Equal variance: Levene's test checks whether variance of total_energy_J
# is similar across the four Treatment x Exercise groups (an assumption
# of standard ANOVA).
leveneTest(total_energy_J ~ Treatment * Exercise,
           data = df_swim)

# Diagnostic plots: set up a 2x2 plotting grid, then draw the standard
# base-R lm/aov diagnostic plots (residuals vs. fitted, Q-Q plot,
# scale-location, residuals vs. leverage) to visually inspect model fit.
par(mfrow = c(2,2))
plot(model_swim)

# -----------------------------
# POST HOC TESTS
# -----------------------------

# emmeans computes the estimated marginal (adjusted) means for every
# combination of Treatment x Exercise, then `pairwise ~` requests all
# pairwise comparisons between those combinations (with multiple-
# comparison correction), to see which specific group pairs differ.
emmeans(model_swim,
        pairwise ~ Treatment * Exercise)

# =========================================================
# WALKING 2-WAY ANOVA
# =========================================================

# Same analysis as above, repeated for the walking subset.
model_walk <- aov(total_energy_J ~ Treatment * Exercise,
                  data = df_walk)

summary(model_walk)

# =========================================================
# WALKING
# =========================================================

perm_walk <- aovperm(
  total_energy_J ~ Treatment * Exercise,
  data = df_walk,
  np = 5000
)

summary(perm_walk)

# -----------------------------
# ASSUMPTION TESTS
# -----------------------------

# Normality
shapiro.test(residuals(model_walk))

# Equal variance
leveneTest(total_energy_J ~ Treatment * Exercise,
           data = df_walk)

# Diagnostic plots
par(mfrow = c(2,2))
plot(model_walk)

# -----------------------------
# POST HOC TESTS
# -----------------------------

emmeans(model_walk,
        pairwise ~ Treatment * Exercise)

# =========================================================
# DEFINE GROUP ORDER
# =========================================================

# `group` is presumably a combined Treatment x Exercise label already
# present in the source CSV (e.g., "TT" = terrestrial+treated,
# "AU" = aquatic+untreated, etc.). Fixing the factor level order here
# controls both the left-to-right order of boxes in the plots below and
# ensures the color mapping (grpColors, next section) lines up correctly
# with the intended group.
group_order <- c("TT", "TU", "AT", "AU")

df_swim$group <- factor(df_swim$group,
                        levels = group_order)

df_walk$group <- factor(df_walk$group,
                        levels = group_order)

# =========================================================
# GROUP COLORS
# =========================================================

# A named vector mapping each group label to a specific hex color, so
# both the swimming and walking plots use a consistent, deliberately
# chosen color scheme (paired dark/light shades for terrestrial vs.
# aquatic-derived groups) rather than ggplot's default palette.
grpColors <- c(
  "TT" = "#8B4513",   # dark brown
  "TU" = "#CC843F",   # tan
  "AT" = "#1A4C99",   # dark blue
  "AU" = "#4DB3FF"    # sky blue
)

# =========================================================
# SWIMMING PLOT
# =========================================================

ggplot(df_swim,
       aes(x = group,
           y = total_energy_J,
           fill = group)) +
  
  # Boxplot summarizing median/IQR per group. outlier.shape = NA hides
  # ggplot's default outlier points so they aren't double-plotted on top
  # of the jittered raw points added next.
  geom_boxplot(outlier.shape = NA,
               alpha = 0.8) +
  
  # Overlay every individual data point, horizontally jittered so
  # overlapping values are visible, giving a sense of sample size and
  # spread on top of the box summary.
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.7) +
  
  # Apply the custom color scheme defined above.
  scale_fill_manual(values = grpColors) +
  
  # Clean, publication-style theme with no gridlines/background panel.
  theme_classic(base_size = 14) +
  
  labs(title = "Swimming Total Energy Excess",
       x = "Group",
       y = "Total Energy Excess (J)") +
  
  # Legend is redundant with the x-axis labels, so it's hidden.
  theme(legend.position = "none")

# =========================================================
# Save plots
# =========================================================

# ggsave() with no explicit `plot =` argument saves the most recently
# displayed ggplot (i.e., the swimming plot built immediately above) as
# a vector SVG file at 5in x 5in.
ggsave(filename = file.path(out_dir,
                            "swimmingEnergyExcessBoxplot.svg"),
       width = 5,
       height = 5
)

# =========================================================
# WALKING PLOT
# =========================================================

# Identical plot construction to the swimming version, but for the
# walking subset.
ggplot(df_walk,
       aes(x = group,
           y = total_energy_J,
           fill = group)) +
  
  geom_boxplot(outlier.shape = NA,
               alpha = 0.8) +
  
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.7) +
  
  scale_fill_manual(values = grpColors) +
  
  theme_classic(base_size = 14) +
  
  labs(title = "Walking Total Energy Excess",
       x = "Group",
       y = "Total Energy Excess (J)") +
  
  theme(legend.position = "none")

# =========================================================
# Save plots
# =========================================================

ggsave(
  filename = file.path(out_dir,
                       "walkingEnergyExcessBoxplot.svg"),
  width = 5,
  height = 5
)

# =========================================================
# SAVE PERMUTATION ANOVA RESULTS
# =========================================================

# Convert the permutation ANOVA summary objects (F-stats, permutation
# p-values, etc.) into plain data frames so they can be written to CSV.
perm_swim_table <- as.data.frame(summary(perm_swim))

perm_walk_table <- as.data.frame(summary(perm_walk))

# row.names = TRUE keeps the term names (Treatment, Exercise,
# Treatment:Exercise) as the first column of the CSV.
write.csv(
  perm_swim_table,
  file.path(out_dir,
            "swimmingEnergyPermANOVA.csv"),
  row.names = TRUE
)

write.csv(
  perm_walk_table,
  file.path(out_dir,
            "walkingEnergyPermANOVA.csv"),
  row.names = TRUE
)

# =========================================================
# POWER CURVE FUNCTION
# PERMUTATION ANOVA (parallelized safely for RStudio/macOS)
# =========================================================
#
# GOAL: estimate statistical power (probability of detecting a true
# effect, if one exists) for the Treatment:Exercise INTERACTION term
# specifically, using the permutation ANOVA, across a range of assumed
# interaction effect sizes. This tells you, e.g., "if the true
# interaction effect were 1.5x as large as what we observed, we'd have
# an estimated XX% chance of detecting it as significant with our
# current sample size."
#
# Approach: parametric bootstrap / simulation-based power analysis --
# fit a model to the REAL data to get realistic coefficients and
# residual noise, then repeatedly simulate NEW synthetic datasets with
# those coefficients scaled up/down, re-run the permutation ANOVA on
# each simulated dataset, and see what fraction come out significant.

library(future)        # backend for parallel/asynchronous evaluation
library(future.apply)  # parallel versions of apply-family functions (future_lapply)

run_power_curve <- function(df,
                            mode_name,
                            nsim = 200,        # number of simulated datasets per effect size
                            np_perm = 299,     # number of permutations per permutation-ANOVA call (kept low here for speed, since this runs nsim x length(effect_sizes) times)
                            n_cores = max(1, parallel::detectCores() - 1)){
  
  # --------------------------------------------------------
  # SET UP PARALLEL WORKERS
  # (PSOCK/multisession -- safe inside RStudio, unlike mclapply)
  # --------------------------------------------------------
  
  # `multisession` spins up separate background R processes to run
  # simulations in parallel. This is used (rather than `mclapply`'s
  # fork-based parallelism) because forking is unreliable/unsafe within
  # RStudio and on macOS.
  plan(multisession, workers = n_cores)
  # Ensure that once this function exits (however it exits, including on
  # error), the parallel workers are shut down and the plan reverts to
  # single-threaded ("sequential") execution, so they don't linger.
  on.exit(plan(sequential), add = TRUE)
  
  # --------------------------------------------------------
  # EFFECT SIZE SCALING
  # --------------------------------------------------------
  
  # Multipliers applied to the *observed* effect (the model coefficients
  # fit to the real data). 1 = same effect size as actually observed;
  # 0.25 = a much weaker hypothetical effect; 3 = a much stronger one.
  # Sweeping across this range traces out how power changes as the true
  # effect gets bigger.
  effect_sizes <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3)
  
  # Accumulator for results (one row per effect size, appended in the loop).
  power_results <- data.frame()
  
  # --------------------------------------------------------
  # BASE MODEL
  # --------------------------------------------------------
  
  # Fit an ordinary linear model to the real data purely to extract:
  #   - coefs: the fitted Treatment/Exercise/interaction effect sizes
  #   - resid_sd: the residual (noise) standard deviation
  # These become the "ground truth" used to generate realistic simulated
  # data below.
  model_lm <- lm(total_energy_J ~ Treatment * Exercise,
                 data = df)
  
  coefs <- coef(model_lm)
  
  resid_sd <- sd(residuals(model_lm))
  
  # --------------------------------------------------------
  # SIMULATION FUNCTION
  # --------------------------------------------------------
  
  # Generates one synthetic dataset: keeps the same design (same
  # Treatment/Exercise group membership as the real data) but replaces
  # the outcome variable with values simulated from the fitted model,
  # after scaling the coefficients by `scale_factor` (the hypothetical
  # effect size) and adding fresh random noise drawn from a normal
  # distribution with the observed residual SD.
  simulate_data <- function(df,
                            coefs,
                            resid_sd,
                            scale_factor){
    
    sim_df <- df
    
    # Build the design/model matrix (intercept + dummy-coded factors +
    # interaction column) corresponding to each row's actual
    # Treatment/Exercise combination.
    X <- model.matrix(~ Treatment * Exercise,
                      data = sim_df)
    
    # Scale every coefficient (intercept and effects alike) by the
    # chosen effect-size multiplier.
    beta_scaled <- coefs * scale_factor
    
    # Predicted mean energy for each row under the scaled effects.
    mu <- X %*% beta_scaled
    
    # Simulated outcome = scaled predicted mean + random normal noise
    # (same noise SD as the real residuals), replacing the real
    # total_energy_J with this synthetic version.
    sim_df$total_energy_J <- as.numeric(mu) +
      rnorm(nrow(sim_df),
            mean = 0,
            sd = resid_sd)
    
    return(sim_df)
  }
  
  # --------------------------------------------------------
  # RUN SIMULATIONS (parallelized across nsim per effect size)
  # --------------------------------------------------------
  
  # Loop over each hypothetical effect-size multiplier.
  for(es in effect_sizes){
    
    cat("\n")
    cat("====================================\n")
    cat(mode_name, "- Effect Size:", es, "\n")
    cat("====================================\n")
    
    # For this effect size, run `nsim` independent simulated-dataset
    # trials IN PARALLEL across the worker processes set up above.
    result_list <- future_lapply(1:nsim, function(i){
      
      # tryCatch guards against any single simulated permutation-ANOVA
      # call failing (e.g., due to a degenerate simulated dataset) --
      # rather than crashing the whole loop, it records the error and
      # moves on, so a handful of failures don't stop the whole power
      # curve from computing.
      tryCatch({
        
        # Generate one synthetic dataset at this effect size.
        sim_df <- simulate_data(df,
                                coefs,
                                resid_sd,
                                scale_factor = es)
        
        # Run the same permutation ANOVA as the real analysis on the
        # simulated data (with a smaller np_perm = 299 permutations,
        # since this needs to happen nsim x length(effect_sizes) times
        # and full precision isn't needed for a power estimate).
        perm_model <- aovperm(
          total_energy_J ~ Treatment * Exercise,
          data = sim_df,
          np = np_perm
        )
        
        # Pull out the permutation p-value specifically for the
        # Treatment:Exercise interaction term -- this is the effect
        # power is being estimated for.
        pval <- summary(perm_model)[
          "Treatment:Exercise",
          "resampled P(>F)"
        ]
        
        # Record whether this simulated trial was "significant" at
        # alpha = 0.05 (1 = yes, 0 = no). No error occurred, so error
        # field is NA.
        list(value = as.numeric(pval < 0.05), error = NA_character_)
        
      }, error = function(e) {
        # If anything above failed, record NA for the result and save
        # the actual error message text (rather than just a failure
        # count) so failures can be diagnosed later.
        list(value = NA_real_, error = conditionMessage(e))
      })
      
    }, future.seed = TRUE,                    # ensures reproducible, statistically sound RNG across parallel workers
    future.packages = c("permuco"))        # explicitly ships the permuco package to each worker process, since they start as blank R sessions
    
    # Extract the significance indicator (0/1/NA) from each of the nsim
    # trial results.
    sig_vec  <- sapply(result_list, function(x) x$value)
    # Extract any error messages that occurred, dropping the NAs (i.e.,
    # keeping only trials that actually failed).
    err_msgs <- sapply(result_list, function(x) x$error)
    err_msgs <- err_msgs[!is.na(err_msgs)]
    
    # Print the ACTUAL underlying error text (not just a count) so we can
    # diagnose instead of guessing.
    if(length(err_msgs) > 0){
      cat("\n---- ", mode_name, " effect_size=", es,
          ": ", length(err_msgs), "/", nsim, " replicates failed ----\n", sep = "")
      # Show the 5 most common distinct error messages and how often
      # each occurred, to spot systematic failure patterns.
      print(head(sort(table(err_msgs), decreasing = TRUE), 5))
    }
    
    # Number of trials that produced a usable (non-NA) result.
    n_valid <- sum(!is.na(sig_vec))
    
    if(n_valid == 0){
      # If every trial failed at this effect size, there's nothing to
      # estimate power from.
      power_est <- NA_real_
      ci <- c(NA_real_, NA_real_)
    } else {
      # Power estimate = proportion of valid trials that were
      # significant at p < .05.
      power_est <- mean(sig_vec, na.rm = TRUE)
      # 95% confidence interval on that proportion, via the exact
      # binomial test (treats "number significant out of n_valid" as a
      # binomial proportion).
      ci <- binom.test(sum(sig_vec, na.rm = TRUE), n_valid)$conf.int
    }
    
    # Append this effect size's results as one row to the results table.
    power_results <- rbind(
      power_results,
      data.frame(
        mode = mode_name,
        effect_size = es,
        power = power_est,
        lower_CI = ci[1],
        upper_CI = ci[2],
        n_valid = sum(!is.na(sig_vec))
      )
    )
  }
  
  return(power_results)
}

# =========================================================
# RUN POWER CURVES (WALKING + SWIMMING)
# =========================================================

# Run the full power-curve simulation (200 simulated datasets x 7 effect
# sizes = 1400 permutation ANOVAs, each with 299 permutations) separately
# for walking and swimming data.
power_walk <- run_power_curve(
  df = df_walk,
  mode_name = "Walking",
  nsim = 200,
  np_perm = 299
)

power_swim <- run_power_curve(
  df = df_swim,
  mode_name = "Swimming",
  nsim = 200,
  np_perm = 299
)

# Combine both modes' results into one long table for plotting/export.
power_all <- rbind(power_walk, power_swim)

# =========================================================
# SAVE POWER CURVE RESULTS
# =========================================================

write.csv(
  power_all,
  file.path(out_dir, "powerCurveInteraction.csv"),
  row.names = FALSE
)

# =========================================================
# PLOT POWER CURVE
# =========================================================

# Line plot of estimated power (y) vs. effect-size multiplier (x),
# one line per Mode, with error bars showing the binomial CI and a
# reference line at the conventional 80% power threshold.
ggplot(power_all,
       aes(x = effect_size,
           y = power,
           color = mode)) +
  
  # Reference line at the conventional "adequately powered" threshold.
  geom_hline(yintercept = 0.8,
             linetype = "dashed",
             color = "grey40") +
  
  geom_line(linewidth = 1) +
  
  geom_point(size = 2) +
  
  # Vertical error bars showing the 95% CI around each power estimate.
  geom_errorbar(aes(ymin = lower_CI,
                    ymax = upper_CI),
                width = 0.05,
                alpha = 0.6) +
  
  # Match line/point colors to the same brown/blue scheme used in the
  # boxplots above, for visual consistency across figures.
  scale_color_manual(values = c("Walking" = "#8B4513",
                                "Swimming" = "#1A4C99")) +
  
  theme_classic(base_size = 14) +
  
  labs(title = "Power Curve: Treatment x Exercise Interaction",
       subtitle = "Dashed line = conventional 0.80 power threshold",
       x = "Effect Size Scaling Factor (1 = observed effect)",
       y = "Power",
       color = "Mode")

ggsave(
  filename = file.path(out_dir, "powerCurveInteractionPo.svg"),
  width = 6,
  height = 5
)