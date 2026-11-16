# RESEARCH QUESTION: Do school conditions causally affect students’ stress levels, and if so, how much and through which pathways?

# ======== I.Install packages =========
install.packages(c(
  "tidyverse",        # data wrangling + plotting
  "janitor",          # tidy column names
  "corrplot",         # correlation heatmap for EDA
  "bnlearn",          # build DAGs, score them with BIC, learn structure
  "dagitty",          # identification: adjustment sets, testable implications
  "ggdag",            # draw DAGs nicely with ggplot
  "MASS",             # polr(): ordinal logistic regression
  "marginaleffects",  # turn model coefficients into interpretable effects
  "sensemakr",        # sensitivity to unmeasured confounding
  "WeightIt",         # optional: inverse-probability weighting (ATE)
  "cobalt"            # optional: check covariate balance after weighting
))

install.packages("BiocManager"); BiocManager::install("Rgraphviz") #for nicer DAG plots

# ---- Load ----
library(tidyverse)
library(janitor)
library(corrplot)
library(bnlearn)
library(dagitty)
library(ggdag)
library(MASS)
library(marginaleffects)
library(sensemakr)
library(WeightIt)
library(cobalt)
select <- dplyr::select

# ========== II.Load & inspect the data ========
# ---- 1. Load the csv and standardize column names ---------------------------------------------------------
stress <- read_csv("/Users/macminhanh/Documents/UChicago/Causal Models/DAG_Student_Stress_Monitoring/Student Stress Monitoring Dataset/StressLevelDataset.csv") |> clean_names()

view(stress)

# ---- 2. Structure: looking at the structure of the data ---------------------------------------------------------
dim(stress)
glimpse(stress)
# Codebook:
# ---- PSYCHOlOGICAL ---------------------------------------------------------
# anxiety_level: Anxiety score. Higher = more anxious
# self_esteem: Self-esteem score. Higher = more self-esteem
# mental_health_history: Prior mental health issue (1 = yes, 0 = no). Binary
# depression: Depression score. Higher = worse
# ---- PHYSIOLOGICAL ---------------------------------------------------------
# headache: Headache frequency. Scale 0-5. Higher = more frequent
# blood_pressure: Blood pressure category ( 1 = low, 2 = normal, 3 = high)
# sleep_quality: Scale 0-5. Higher = better sleep
# breathing_problem: Scale 0-5. Breathing difficulty. Higher = more difficult
# --- ENVIRONMENTAL ---------------------------------------------------
# noise_level: Scale 0-5. Noise where the student lives/studies. Higher = noisier
# living_conditions: Scale 0-5. Higher = better
# safety: Scale 0-5. Feeling safe. Higher = better
# basic_needs: Scale 0-5. Basic needs met (food, shelter). Higher = better
# --- ACADEMIC --------------------------------------------------------
# academic_performance: Scale 0-5.  Higher = better grades
# study_load: Scale 0-5.  Higher = heavier workload
# teacher_student_relationship: Scale 0-5.  Higher = better relationship
# future_career_concerns: Scale 0-5.  Higher = more worried
# --- SOCIAL ----------------------------------------------------------
# social_support: Higher = more support (check range)
# peer_pressure: Scale 0-5. Higher = more pressure
# extracurricular_activities: Scale 0-5. Level of involvement
# bullying: Scale 0-5. Higher = more bullying
#--- OUTCOME ---------------------------------------------------------
# stress_level: 0 = low, 1 = medium, 2 = high (ordinal)

# ---- 3. Missing values ---------------------------------------------------------
colSums(is.na(stress)) # no missing value

# ---- 4. Range of every varaible ---------------------------------------------------------
stress |>
  summarise(across(everything(), list(min = min, max = max))) |>
  pivot_longer(everything(),
               names_to = c("variable", ".value"),
               names_pattern = "(.*)_(min|max)") |> 
  print(n = Inf)

# ---- 5. Check if the sample is imbalance ---------------------------------------------------------
stress |> count(stress_level) |> mutate(pct = round(100 * n / sum(n), 1))
#  0   373  33.9
#  1   358  32.5
#  2   369  33.5

# ---- 6. Distribution of every variable in one picture ---------------------------------------------------------
stress |>
  pivot_longer(everything()) |>
  ggplot(aes(value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  facet_wrap(~ name, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of each variable")

# ---- 7. Correlation heatmap, ordered so related variables cluster together ---------------------------------------------------------
cor_mat <- cor(stress, method = "spearman")   # Spearman: good for ordinal scales
corrplot(cor_mat, method = "color", type = "upper", order = "hclust",
         tl.cex = 0.5, tl.col = "black")

# ---- 8. Which variables are most correlated with stress? ---------------------------------------------------------
cor_mat[, "stress_level"] |> sort() |> round(2)

# ======== III.DAGs =========
vars_keep <- c("mental_health_history", "basic_needs",
               "study_load", "teacher_student_relationship", "bullying",
               "peer_pressure", "social_support",
               "academic_performance", "future_career_concerns",
               "self_esteem", "anxiety_level", "depression", "sleep_quality",
               "stress_level")
df_s <- df[, vars_keep]

# ---- 1.Tiers (earliest -> latest). Arrows may only go forward ---------------------------------------------------------
tiers_s <- list(
  c("mental_health_history", "basic_needs"),                        # T1 Background
  c("study_load", "teacher_student_relationship",
    "bullying", "peer_pressure", "social_support"),                # T2 School & social conditions
  c("academic_performance", "future_career_concerns"),             # T3 Academic outcomes
  c("self_esteem", "anxiety_level", "depression", "sleep_quality"), # T4 Psychological
  c("stress_level")                                                # T5 Outcome
)

# Sanity check: every column is in exactly one tier 
setdiff(names(df_s), unlist(tiers_s))   # should be character(0)

# Forbidden arrows = every backward arrow
bl_s <- tiers2blacklist(tiers_s)
nrow(bl_s)

#---------- DAG 1 ---------------------------------------------------------
#Hypothesis: workload causes stress through sleep quality. I compare DAG A (workload → sleep → stress) with DAG B (which adds a direct workload → stress arrow). 
#Each DAG’s BIC is the sum of its node scores: each variable is regressed on its parents, its log-likelihood measures fit, and a penalty is subtracted for each parameter. 
#The DAG with the higher BIC explains the data better relative to its complexity. BIC can’t determine arrow direction, which comes from theory.

set.seed(42)
dag1 <- hc(df_s, score = "bic-g", blacklist = bl_s)

narcs(dag1)                              # how many arrows
bnlearn::parents(dag1, "stress_level")   # what points directly into stress
score(dag1, df, type = "bic-g")          # its BIC (higher = better in bnlearn)

# Plot
bn_to_dagitty <- function(bn) {
  a <- arcs(bn)
  dagitty(paste0("dag {", paste(a[, 1], "->", a[, 2], collapse = "; "), "}"))
}

ggdag(bn_to_dagitty(dag1), layout = "sugiyama", text = FALSE, use_labels = "name") +
  theme_dag()

stress_parents <- bnlearn::parents(dag1, "stress_level")
focus <- dagitty(paste0("dag {", paste(stress_parents, "-> stress_level", collapse = "; "), "}"))
ggdag(focus, layout = "star", text = FALSE, use_labels = "name") + theme_dag()

# Very dense arrow pattern (67 arrows).
# There are 9 direct paretns of stress_level: self-esteem, depression, sleep quality, basic needs, academic performance, study load, peer pressure, extracurriculars, bullying.
# Most of the variabels point directly to stress
# BIC Score: -24748.98

#---------- DAG 2 ---------------------------------------------------------
# Hypothesis: School conditions affect stress only through students’ psychological state.
# Workload, teacher relationships and career worries don’t raise stress directly. They raise it because they worsen sleep, anxiety, depression or self-esteem.

# Extra rule: nothing from tiers 1-3 may point directly at stress_level
upstream_s <- unlist(tiers_s[1:3])
bl_med_s <- rbind(as.data.frame(bl_s),
                  data.frame(from = upstream_s, to = "stress_level"))

dag2 <- hc(df_s, score = "bic-g", blacklist = bl_med_s)

narcs(dag2)                              # number of arrows
bnlearn::parents(dag2, "stress_level")   # should be psych variables only
score(dag2, df_s, type = "bic-g")        # BIC of DAG 2
score(dag1, df_s, type = "bic-g")        # BIC of DAG 1, for comparison

ggdag(bn_to_dagitty(dag2), layout = "sugiyama", text = FALSE, use_labels = "name") +
  theme_dag()

stress_parents <- bnlearn::parents(dag2, "stress_level")
focus <- dagitty(paste0("dag {", paste(stress_parents, "-> stress_level", collapse = "; "), "}"))
ggdag(focus, layout = "star", text = FALSE, use_labels = "name") + theme_dag()

# Dense arrow pattern (63 arrows).
# There are only 4 direct paretns of stress_level: self-esteem, depression, sleep quality, and anxiety_level
# BIC Score: -24845.25.

#---------- DAG 3 ---------------------------------------------------------
#Hypothesis: Stress comes first, and anxiety, depression, poor sleep and low self-esteem are results of stress, not causes.
# School conditions create stress directly, and that stress then shows up as anxiety, low mood, bad sleep and lower self-esteem.

tiers_rev_s <- list(
  tiers_s[[1]],      # T1 Background
  tiers_s[[2]],      # T2 School & social conditions
  tiers_s[[3]],      # T3 Academic outcomes
  "stress_level",    # T4 Outcome (moved up)
  tiers_s[[4]]       # T5 Psych (now AFTER stress)
)
bl_rev_s <- tiers2blacklist(tiers_rev_s)

dag3 <- hc(df_s, score = "bic-g", blacklist = bl_rev_s)

narcs(dag3)
bnlearn::parents(dag3, "stress_level")    # should contain no psych variables
bnlearn::children(dag3, "stress_level")   # psych variables stress points to
score(dag3, df_s, type = "bic-g")

ggdag(bn_to_dagitty(dag3), layout = "sugiyama", text = FALSE, use_labels = "name") +
  theme_dag()

stress_parents <- bnlearn::parents(dag3, "stress_level")
focus <- dagitty(paste0("dag {", paste(stress_parents, "-> stress_level", collapse = "; "), "}"))
ggdag(focus, layout = "star", text = FALSE, use_labels = "name") + theme_dag()

# Dense arrow pattern (63 arrows).
# There are 7 direct paretns of stress_level but none of them are psych variables
# BIC Score: -24737.66

#---------- DAG 4 ---------------------------------------------------------
# Hypothesis: Poor sleep and low self-esteem cause stress, but anxiety and depression are results of stress.
# Sleep and self-esteem are fairly stable, everyday conditions that make a student more vulnerable to stress. Anxiety and depression are emotional reactions that build up after a student is stressed.

tiers_mix_s <- list(
  tiers_s[[1]],                         # T1 Background
  tiers_s[[2]],                         # T2 School & social conditions
  tiers_s[[3]],                         # T3 Academic outcomes
  c("sleep_quality", "self_esteem"),    # T4 Psych causes
  "stress_level",                       # T5 Outcome
  c("anxiety_level", "depression")      # T6 Psych results
)
setdiff(names(df_s), unlist(tiers_mix_s))   # should be character(0)
bl_mix_s <- tiers2blacklist(tiers_mix_s)

dag4 <- hc(df_s, score = "bic-g", blacklist = bl_mix_s)

narcs(dag4)
bnlearn::parents(dag4, "stress_level")    # no anxiety/depression here
bnlearn::children(dag4, "stress_level")   # anxiety/depression may appear here
score(dag4, df_s, type = "bic-g")

# Dense arrow pattern (68 arrows).
# There are 7 direct paretns of stress_level but none of them are anxiety or depression related
# BIC Score: -24745.37

# Compare 4 dag models
dags_s <- list(learned = dag1, mediated = dag2, reverse = dag3, mixed = dag4)

comparison <- tibble(
  dag   = names(dags_s),
  arcs  = sapply(dags_s, narcs),
  bic_g = sapply(dags_s, function(d) score(d, df_s, type = "bic-g")),
  stress_parents = sapply(dags_s, function(d)
    paste(bnlearn::parents(d, "stress_level"), collapse = ", "))
) |>
  arrange(desc(bic_g)) |>
  mutate(gap_from_best = bic_g - max(bic_g))

comparison 
comparison|> select(dag, stress_parents) |> print(width = Inf)

# Hypothesis 3(Stress comes first, and anxiety, depression, poor sleep and low self-esteem are results of stress, not causes) has the highest BIC.

# ------- Rerun with random restarts -----------------------------
#hc() starts with an empty diagram and keeps adding, removing or flipping one arrow at a time, whichever raises the BIC most, until nothing helps, while never breaking the rules.
# Re-learn each DAG with 20 random restarts to avoid local optima
set.seed(42)
dag1_r <- hc(df_s, score = "bic-g", blacklist = bl_s,     restart = 20, perturb = 5)
dag2_r <- hc(df_s, score = "bic-g", blacklist = bl_med_s, restart = 20, perturb = 5)
dag3_r <- hc(df_s, score = "bic-g", blacklist = bl_rev_s, restart = 20, perturb = 5)
dag4_r <- hc(df_s, score = "bic-g", blacklist = bl_mix_s, restart = 20, perturb = 5)

dags_r <- list(learned = dag1_r, mediated = dag2_r, reverse = dag3_r, mixed = dag4_r)
tibble(dag   = names(dags_r),
       bic_g = sapply(dags_r, function(d) score(d, df_s, type = "bic-g"))) |>
  arrange(desc(bic_g)) |>
  mutate(gap_from_best = bic_g - max(bic_g))

#Dag 3 is still the best one.

# ------- Let R draw everything itself -----------------------------
set.seed(42)
dag_free <- hc(df_s, score = "bic-g", restart = 20, perturb = 5)

narcs(dag_free)
score(dag_free, df_s, type = "bic-g")              # compare with Reverse
bnlearn::parents(dag_free, "stress_level")
bnlearn::children(dag_free, "stress_level")

# Arrows R drew that break our tiers (i.e. "nonsense" by our theory)
violations <- merge(as.data.frame(arcs(dag_free)), as.data.frame(bl_s))
nrow(violations)
violations

#BIC score: -24720.18 -> 
# The unconstrained DAG contained 37 arcs violating the theorized ordering. Twelve pointed into background variables (e.g., anxiety → mental health history), which is temporally impossible, confirming that BIC alone cannot establish causal direction. Notably, the unconstrained search independently oriented stress → anxiety, depression and sleep quality, consistent with the best-fitting Reverse DAG. Some violations (e.g., stress → academic performance; psychological states → perceived school conditions) suggest plausible reverse causation or self-report perception bias, which we note as limitations.

# ---- Mixed-type BIC----------------
# bic-g treated every variable as continuous. The custom BIC scores the same four DAGs, but uses the model that fits each variable’s real type.
mixed_bic <- function(dag, data,
                      binary  = "mental_health_history",
                      ordinal = "stress_level") {
  node_scores <- sapply(bnlearn::nodes(dag), function(v) {
    pa <- bnlearn::parents(dag, v)
    f  <- reformulate(if (length(pa)) pa else "1", response = v)
    
    fit <- if (v %in% binary) {
      glm(f, data = data, family = binomial)                   # logistic
    } else if (v %in% ordinal) {
      polr(update(f, factor(.) ~ .), data = data, Hess = TRUE) # ordinal
    } else {
      lm(f, data = data)                                       # linear
    }
    -BIC(fit) / 2      # convert to bnlearn's scale: higher = better
  })
  sum(node_scores)
}

# Score the four (restart) DAGs with both judges
tibble(
  dag       = names(dags_r),
  bic_g     = sapply(dags_r, function(d) score(d, df_s, type = "bic-g")),
  bic_mixed = sapply(dags_r, function(d) mixed_bic(d, df_s))
) |>
  mutate(gap_g     = bic_g - max(bic_g),
         gap_mixed = bic_mixed - max(bic_mixed)) |>
  arrange(desc(bic_mixed))

# Both the Gaussian BIC and a mixed-type BIC (logistic for binary, ordinal logistic for stress level, linear otherwise) ranked the DAGs identically: Reverse > Mixed > Learned > Mediated. The Reverse DAG, in which stress precedes anxiety, depression, sleep quality and self-esteem, was selected as the final causal model.

# ======== IV.Identification with the Reverse DAG =========
g_rev <- bn_to_dagitty(dag3_r)

causes <- c("study_load", "teacher_student_relationship", "future_career_concerns")

# --- 1. Is there any causal path from each cause to stress? ---------
for (x in causes) {
  cat("\n====", x, "====\n")
  p <- paths(g_rev, from = x, to = "stress_level", directed = TRUE)
  cat("Directed (causal) paths:", sum(p$open), "\n")
  print(head(p$paths, 5))
}

# All three school factors have causal paths to stress in the Reverse DAG, so all three effects are worth estimating
# Each school factor had both a direct path and indirect paths to stress, mainly through future career concerns and academic performance. The teacher–student relationship operated largely through school-climate variables (bullying, peer pressure); because these share a tier, their relative orientation was not determined by theory, and we assess sensitivity to it.

"teacher_student_relationship" %in% bnlearn::parents(dag3_r, "stress_level")

# ---- 2. Adjustment sets for the total effect (backdoor criterion) ------
for (x in causes) {
  cat("\n==== Total effect of", x, "====\n")
  print(adjustmentSets(g_rev, exposure = x, outcome = "stress_level",
                       effect = "total"))
}

# ---- 3. Adjustment sets for the DIRECT effect (not via academic mediators)------
for (x in causes) {
  cat("\n==== Direct effect of", x, "====\n")
  print(adjustmentSets(g_rev, exposure = x, outcome = "stress_level",
                       effect = "direct"))
} 

# ========= V. Estimation ==================================
# Naive effect: with no controls at all -> misleading. For comparison only
# Total effect: If a school changed this condition, how much would stress chaneg overall? Counts every path from the cause to stress. Add confounders
# Direct effect: How much of the effect goes straight to stress, without passing through other school factors? Adding mediators here
# Indirect effect: How much of the effect travels through other variables (the mediators)? Indirect effect = total - direct effect

# ---- 1. Total, direct, indirect effects -------------------
# Adjustment sets from Step 4 (Reverse DAG)
adj_total <- list(
  study_load = c("bullying", "mental_health_history", "peer_pressure", "social_support"),
  teacher_student_relationship = c("basic_needs", "mental_health_history"),
  future_career_concerns = c("basic_needs", "bullying", "mental_health_history",
                             "peer_pressure", "social_support", "study_load",
                             "teacher_student_relationship")
)

adj_direct <- list(
  study_load = c("academic_performance", "basic_needs", "bullying",
                 "future_career_concerns", "mental_health_history",
                 "peer_pressure", "social_support"),
  teacher_student_relationship = c("academic_performance", "basic_needs", "bullying",
                                   "future_career_concerns", "mental_health_history",
                                   "peer_pressure", "social_support", "study_load"),
  future_career_concerns = c("academic_performance", "basic_needs", "bullying",
                             "mental_health_history", "peer_pressure",
                             "social_support", "study_load")
)

# Standardized copy (mean 0, SD 1) so the three causes can be compared
df_z <- df_s |> mutate(across(everything(), ~ as.numeric(scale(.))))

# One regression: cause + controls -> estimate, 95% CI, p-value
fit_effect <- function(x, controls, data) {
  m  <- lm(reformulate(c(x, controls), "stress_level"), data = data)
  ci <- confint(m)[x, ]
  c(est = unname(coef(m)[x]), lo = unname(ci[1]), hi = unname(ci[2]),
    p = summary(m)$coefficients[x, 4])
}

# Table of naive, total, direct and indirect effects
effects_table <- function(data) {
  map_dfr(names(adj_total), function(x) {
    naive <- coef(lm(reformulate(x, "stress_level"), data = data))[x]
    tot   <- fit_effect(x, adj_total[[x]],  data)
    dir   <- fit_effect(x, adj_direct[[x]], data)
    tibble(cause     = x,
           naive     = unname(naive),
           total     = tot["est"], total_lo  = tot["lo"], total_hi  = tot["hi"],
           direct    = dir["est"], direct_lo = dir["lo"], direct_hi = dir["hi"],
           indirect  = tot["est"] - dir["est"])
  })
}

effects_raw <- effects_table(df_s)   # raw units
effects_std <- effects_table(df_z)   # standardized

effects_raw |> mutate(across(where(is.numeric), ~ round(.x, 3))) |> print(width = Inf)
effects_std |> mutate(across(where(is.numeric), ~ round(.x, 3))) |> print(width = Inf)

# Robustness: study_load direct effect using the second valid adjustment set
fit_effect("study_load",
           c("basic_needs", "bullying", "future_career_concerns", "mental_health_history",
             "peer_pressure", "social_support", "teacher_student_relationship"),
           df_s) |> round(3)

# RESULTS: All three school factors have significant total effects on stress (p < .001),
# in the expected direction. Raw: study load +0.10 [0.08, 0.13], career concerns +0.09
# [0.06, 0.12], teacher relationship -0.17 [-0.20, -0.14] stress levels per 1-point increase.
# Standardized: teacher relationship has the largest effect (beta = -0.29), vs career
# concerns (0.17) and study load (0.16). Naive estimates were ~2-4x larger, so most of the
# raw correlation was confounding (mental health history, basic needs, social climate).
# Study load and career concerns act mostly directly; the teacher relationship has NO direct
# effect (-0.003, CI includes 0, as the DAG predicted). Its effect works entirely through
# school climate (bullying, peer pressure, career concerns, grades).
# Caveat: the teacher effect depends on the teacher -> bullying direction (tested in 5c).

# ---- 2. Ordinal logistic regression (total effects) ----------------
# Stress is not a real number, it is three ordered categories (low, medium, and high)
# lm assumes equal spacing, treatign the jump fro low to medium and medium to high the same size. We don't know how far the distance between each categories in real life
# Ordinal logistic regression predicts the probabilities for each stress level
# Each point of (ex: better teacher relationship) lowers the chance of high stress by 10 percentage points.
# We use both lm and polr since lm is easy to compute total, direct and indirect effects. We use polr for roburthness check.

ordinal_effects <- map_dfr(names(adj_total), function(x) {
  d <- df_s |> mutate(stress_f = factor(stress_level, ordered = TRUE))
  m <- polr(reformulate(c(x, adj_total[[x]]), "stress_f"), data = d, Hess = TRUE)
  
  b  <- coef(m)[x]
  se <- sqrt(diag(vcov(m)))[x]
  s  <- as.data.frame(avg_slopes(m, variables = x))   # effect on each category's probability
  
  tibble(cause          = x,
         log_odds       = unname(b),
         odds_ratio     = unname(exp(b)),
         z              = unname(b / se),
         p_value        = 2 * pnorm(-abs(unname(b / se))),
         dP_low_stress  = s$estimate[s$group == "0"],
         dP_high_stress = s$estimate[s$group == "2"])
})

ordinal_effects |>
  mutate(across(c(log_odds, odds_ratio, z, dP_low_stress, dP_high_stress), ~ round(.x, 3)),
         p_value = format.pval(p_value, digits = 3, eps = 0.001))

# RESULTS (ordinal logistic, polr): All three total effects keep the same sign and
# remain significant (p < .001), so section 1 does not depend on treating stress as continuous.
# Odds ratios: study load 1.77, career concerns 1.64, teacher relationship 0.42 (each
# 1-point improvement cuts the odds of a higher stress category by 58%). Per 1-point
# increase, P(high stress) changes by +3.8 pp (study load), +3.3 pp (career concerns)
# and -10.1 pp (teacher relationship). The teacher relationship remains the strongest
# effect, matching the section 1 ranking.

# ----- 3. DAG sensitivity check --------------------------
# DAG sensitivity check is needed because every causal estimate depends on my DAG hypothesis. IF my DAG hypothesis is wrong, the adjustment sets are wrong, so are the effects.
# My Reverse DAG won but the Mixed DAG was only abotu 7 BIC points behind. 
# Using hc(), it chooses teacher -> bullying but the data scored the opposite direction the same. The data can't decide which direction is right.
# Here, DAG rebuild the path under each plausible alternative, let dagity find the new adjustment sets automatically, re-estimate each effect and then compare

# (1) Flipped DAG: Reverse, but bullying & peer pressure come BEFORE teacher relationship
bl_flip <- rbind(as.data.frame(bl_rev_s),
                 data.frame(from = "teacher_student_relationship",
                            to   = c("bullying", "peer_pressure")))
set.seed(42)
dag3_flip <- hc(df_s, score = "bic-g", blacklist = bl_flip, restart = 20, perturb = 5)

# Does the data prefer one direction? (small gap = data can't tell)
score(dag3_r,    df_s, type = "bic-g")
score(dag3_flip, df_s, type = "bic-g")
# Each one-point improvement in the teacher–student relationship reduces stress by between 0.09 and 0.17 levels, depending on whether bullying and peer pressure are assumed to follow or precede the relationship. The two orderings fit the data equally well (ΔBIC = 0.16).

# Helper: for any DAG, find adjustment sets and estimate total effects
effects_for_dag <- function(dag, label) {
  g <- bn_to_dagitty(dag)
  map_dfr(causes, function(x) {
    sets <- adjustmentSets(g, exposure = x, outcome = "stress_level", effect = "total")
    if (length(sets) == 0) {
      return(tibble(dag = label, cause = x, adj_set = "NOT IDENTIFIED",
                    est = NA_real_, lo = NA_real_, hi = NA_real_))
    }
    z <- as.character(sets[[1]])
    r <- fit_effect(x, z, df_s)
    tibble(dag = label, cause = x, adj_set = paste(z, collapse = ", "),
           est = r["est"], lo = r["lo"], hi = r["hi"])
  })
}

# (2) Compare: Reverse (main), Flipped, Mixed, Learned
sens <- bind_rows(
  effects_for_dag(dag3_r,    "reverse (main)"),
  effects_for_dag(dag3_flip, "reverse + bullying->teacher"),
  effects_for_dag(dag4_r,    "mixed"),
  effects_for_dag(dag1_r,    "learned")
)

# Side-by-side estimates
sens |>
  select(dag, cause, est) |>
  pivot_wider(names_from = dag, values_from = est) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Full detail (adjustment sets + CIs)
sens |> mutate(across(c(est, lo, hi), ~ round(.x, 3))) |> print(n = Inf, width = Inf)

# RESULTS (DAG sensitivity): Study load (+0.101) and career concerns (+0.089) are
# identical under all four DAGs (Reverse, Reverse + bullying->teacher, Mixed, Learned):
# fully robust. The teacher relationship effect depends on the DAG: -0.173 (Reverse),
# -0.086 (bullying & peer pressure before teacher), -0.104 (Mixed, Learned). When bullying
# is treated as a confounder, the effect roughly halves, but it remains negative and
# significant in every DAG (CIs exclude 0). Conclusion: teacher relationships reduce stress
# Reverse vs Reverse + bullying->teacher fit equally well (BIC -24737.45 vs -24737.61,
# gap = 0.16), so the data cannot determine the direction. Report the teacher effect
# as a range: -0.09 to -0.17 stress levels per 1-point improvement (both significant).
# Study load (+0.10) and career concerns (+0.09) are unaffected.

#---- 4. Sensitivity to unmeasured confounding ---------------
# The backdoor adjustment only works if I have controlled for every confounders. I can only control for variables in my dataset. Something unmeasured could be affecting both a school factor and stress (ex: personality, general health,etc)
# Question asked: How strong would a hidden confounders need to be to wipe out this effect
# Benchmark: What if the hidden confounder were 1x, 2x, 3x as strong as mental health history?

# Adjustment set for teacher under the flipped DAG (from 5c)
adj_teacher_flip <- sens |>
  filter(dag == "reverse + bullying->teacher",
         cause == "teacher_student_relationship") |>
  pull(adj_set) |>
  strsplit(", ") |> unlist()
adj_teacher_flip   # check it

run_sense <- function(x, controls, label) {
  m <- lm(reformulate(c(x, controls), "stress_level"), data = df_s)
  sensemakr(model = m, treatment = x,
            benchmark_covariates = "mental_health_history", kd = 1:3)
}

sens_list <- list(
  "study_load"           = run_sense("study_load", adj_total$study_load),
  "career_concerns"      = run_sense("future_career_concerns", adj_total$future_career_concerns),
  "teacher (main)"       = run_sense("teacher_student_relationship", adj_total$teacher_student_relationship),
  "teacher (flipped)"    = run_sense("teacher_student_relationship", adj_teacher_flip)
)

# Summary: estimate, partial R2, robustness values
imap_dfr(sens_list, function(s, lab) {
  st <- s$sensitivity_stats
  tibble(model      = lab,
         estimate   = st$estimate,
         partial_R2 = st$r2yd.x,     # share of stress variance uniquely explained by the cause
         RV         = st$rv_q,       # strength needed to erase the effect
         RV_sig     = st$rv_qa)      # strength needed to make it non-significant
}) |> mutate(across(where(is.numeric), ~ round(.x, 3)))

# Benchmarks: adjusted estimate if a hidden confounder were 1x, 2x, 3x mental_health_history
imap_dfr(sens_list, function(s, lab) {
  s$bounds |> as_tibble() |> mutate(model = lab, .before = 1)
}) |>
  select(model, bound_label, adjusted_estimate, adjusted_lower_CI, adjusted_upper_CI) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  print(n = Inf)

# RESULTS (sensemakr, benchmark = mental_health_history): Study load (0.101 -> 0.047)
# and career concerns (0.089 -> 0.067) stay positive and significant even with a hidden
# confounder 3x as strong as mental health history: robust. Teacher relationship (main
# DAG, -0.173) is fragile: a confounder 2x as strong makes it non-significant. The flipped
# DAG estimate (-0.086, already adjusting for bullying & peer pressure) survives 3x
# (-0.045, CI excludes 0). Conclusion: the conservative teacher effect (~ -0.09) is robust;
# the larger estimate may be inflated by confounding.

# RV: a hidden confounder would need to explain (of the residual variance in both the
# cause and stress) 20.3% to erase study load (15.4% to make it non-sig), 17.4% for career
# concerns (12.3%), 28.9% for teacher-main (24.6%), 16.3% for teacher-flipped (11.2%).
# Teacher-main has the highest RV, but the benchmark shows confounders 2x as strong as
# mental health history are plausible for it, so it is fragile. Study load, career concerns
# and teacher-flipped survive both checks.

# Ordinal model: Does the result depend on how stress is measured?
# DAG sensitivity: Does the result depend on the DAG we chose?
# Unmeasured confounding: Does the result depend on variables we didn't measure? 

