# =============================================================================
# Polypterus body size: two-way ANOVAs (Environment x Ramp) on mass and length,
# run separately for each experiment x test so every fish appears only once.
#
# Decision rule for each response:
#   1. Fit on raw data; check Shapiro-Wilk on residuals and Levene across groups.
#   2. If either fails (p < 0.05), refit on ln-transformed data.
#   3. If ln data still have unequal variances, use HC3-robust (White-adjusted)
#      ANOVA on the ln data.
# All ANOVAs use Type II sums of squares (design is unbalanced).
# =============================================================================

# install.packages(c("readxl", "car"))   # run once if needed
library(readxl)
library(car)

# ---- Data -------------------------------------------------------------------
file <- "polypterus_all_experiments_2.xlsx"   # change path if needed
d <- read_excel(file, sheet = "All experiments")
d <- as.data.frame(d)
names(d)[names(d) == "Experiment type"] <- "Experiment"
names(d)[names(d) == "Fish ID"]         <- "Fish"
names(d)[names(d) == "Mass (g)"]        <- "Mass"
names(d)[names(d) == "Length (cm)"]     <- "Length"
d <- d[!is.na(d$Experiment), ]

d$Environment <- factor(d$Environment)
d$Ramp        <- factor(d$Ramp, levels = c("NR", "R"))   # NR = untrained, R = trained
d$lnMass      <- log(d$Mass)
d$lnLength    <- log(d$Length)

alpha <- 0.05

# ---- Helper: fit one model and return checks + ANOVA table -------------------
fit_anova <- function(dat, y, robust = FALSE) {
  f   <- as.formula(paste(y, "~ Environment * Ramp"))
  mod <- lm(f, data = dat)
  tab <- if (robust) Anova(mod, type = 2, white.adjust = "hc3")
         else        Anova(mod, type = 2)
  sw  <- shapiro.test(residuals(mod))$p.value
  lev <- leveneTest(f, data = dat)[1, "Pr(>F)"]   # median-centred (Brown-Forsythe)
  list(model = mod, table = tab, shapiro_p = sw, levene_p = lev)
}

# ---- Run the decision rule for one response in one test ---------------------
analyse <- function(dat, var) {
  steps <- list()
  raw <- fit_anova(dat, var)
  steps[["raw"]] <- raw
  final <- raw; label <- paste(var, "(raw, OLS)")

  if (raw$shapiro_p < alpha || raw$levene_p < alpha) {
    lnvar <- paste0("ln", var)
    lg <- fit_anova(dat, lnvar)
    steps[["ln"]] <- lg
    final <- lg; label <- paste0("ln(", var, ") (OLS)")

    if (lg$levene_p < alpha || lg$shapiro_p < alpha) {
      rb <- fit_anova(dat, lnvar, robust = TRUE)
      steps[["ln_hc3"]] <- rb
      final <- rb; label <- paste0("ln(", var, ") (HC3-robust)")
    }
  }
  list(steps = steps, final = final, label = label)
}

# ---- Loop over experiment x test ---------------------------------------------
tests <- list(c("Endurance", "Swim"),
              c("Endurance", "Walk"),
              c("Kinematics", "Swim"),
              c("Kinematics", "Walk"),
              c("Kinematics", "Full terrestrial"))

summary_rows <- list()

for (tt in tests) {
  sub <- d[d$Experiment == tt[1] & d$Test == tt[2], ]
  stopifnot(!any(duplicated(sub$Fish)))   # each fish once per analysis

  cat("\n=====================================================\n")
  cat(tt[1], "-", tt[2], "| n =", nrow(sub), "fish\n")
  cat("=====================================================\n")
  print(table(sub$Environment, sub$Ramp))

  for (v in c("Mass", "Length")) {
    res <- analyse(sub, v)

    cat("\n---", v, "---\n")
    for (s in names(res$steps)) {
      st <- res$steps[[s]]
      cat(sprintf("[%s] Shapiro p = %.3f | Levene p = %.3f\n",
                  s, st$shapiro_p, st$levene_p))
    }
    cat("Model reported:", res$label, "\n")
    print(res$final$table)

    tab <- res$final$table
    for (term in c("Environment", "Ramp", "Environment:Ramp")) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Experiment = tt[1], Test = tt[2], n = nrow(sub),
        Response = v, Model = res$label, Term = term,
        df = tab[term, "Df"], df_error = tab["Residuals", "Df"],
        F = round(tab[term, intersect(c("F value", "F"), colnames(tab))[1]], 2), p = round(tab[term, "Pr(>F)"], 4),
        Shapiro_p = round(res$final$shapiro_p, 3),
        Levene_p  = round(res$final$levene_p, 3)
      )
    }
  }
}

results <- do.call(rbind, summary_rows)
results$sig <- ifelse(results$p < alpha, "*", "")

cat("\n\n================ SUMMARY (models to report) ================\n")
print(results, row.names = FALSE)

# ---- Group means and SDs by test --------------------------------------------
means <- aggregate(cbind(Mass, Length) ~ Experiment + Test + Environment + Ramp,
                   data = d, FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
means <- do.call(data.frame, means)
means <- means[order(means$Experiment, means$Test, means$Environment, means$Ramp), ]
cat("\n================ GROUP MEANS ================\n")
print(means, digits = 3, row.names = FALSE)

# ---- Save -------------------------------------------------------------------
write.csv(results, "size_anova_results.csv", row.names = FALSE)
write.csv(means,   "size_group_means.csv",   row.names = FALSE)
