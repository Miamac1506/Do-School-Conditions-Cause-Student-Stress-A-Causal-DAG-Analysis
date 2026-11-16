# Do School Conditions Cause Student Stress? A Causal DAG Analysis

**Research question:** Do school conditions causally affect students' stress levels? If so, how much, and through which pathways?

This project uses **causal directed acyclic graphs (DAGs)** to go beyond correlation. It learns and compares competing causal structures from the data, finds the right variables to control for using the **backdoor criterion**, estimates total, direct, and indirect effects, and checks how sensitive the results are to the model's assumptions and to unmeasured confounding.

---

## Key Findings

| School factor (1-point increase) | Total effect on stress | Direct effect | Robust? |
|---|---|---|---|
| **Study load** | **+0.10** [0.08, 0.13] | Mostly direct | ✅ Same across all DAGs and robust to hidden confounding |
| **Future career concerns** | **+0.09** [0.06, 0.12] | Mostly direct | ✅ Same across all DAGs and robust to hidden confounding |
| **Teacher–student relationship** | **−0.09 to −0.17** | None (−0.003, CI includes 0) | ⚠️ Size depends on the DAG; the conservative estimate (−0.09) is robust |

*Effects are in stress levels on a 0–2 scale, estimated with linear regression on the total-effect adjustment sets.*

- **All three school conditions cause stress, in the expected direction.** Heavier workloads and career worries raise stress, and better teacher relationships lower it.
- **The teacher–student relationship has the largest effect** (standardized β = −0.29, vs. 0.17 for career concerns and 0.16 for study load). In the ordinal model, each 1-point improvement **lowers the chance of high stress by about 10 percentage points**.
- **The teacher effect is entirely indirect.** It works through the school climate (bullying, peer pressure), career concerns, and grades, not directly on stress.
- **Naive correlations overstate the effects 2–4×.** Most of the raw association comes from confounders such as mental health history, basic needs, and the social climate.
- **Stress comes before psychological symptoms.** The best-fitting DAG has stress causing anxiety, depression, poor sleep, and low self-esteem, not the other way around.

---

## Dataset

The **Student Stress Monitoring Dataset** (`StressLevelDataset.csv`) has 1,100 students, 21 survey measures, and no missing values.

- **Outcome:** `stress_level`, an ordered category: 0 = low, 1 = medium, 2 = high. The classes are balanced (34% / 33% / 34%).
- **Predictors** fall into five groups: psychological, physiological, environmental, academic, and social.

13 of these variables are used in the causal model, arranged in **time-ordered tiers**. Arrows may only point forward in time:

| Tier | Variables |
|---|---|
| T1 Background | mental health history, basic needs |
| T2 School & social conditions | study load, teacher–student relationship, bullying, peer pressure, social support |
| T3 Academic outcomes | academic performance, future career concerns |
| T4 Psychological | self-esteem, anxiety, depression, sleep quality |
| T5 Outcome | stress level |

---

## Methodology

### 1. EDA
Variable distributions, a Spearman correlation heatmap (suited to ordinal scales), and each variable's correlation with stress.

### 2. Learning and Comparing Causal DAGs
Four competing hypotheses were encoded as **blacklists** (forbidden arrows) and learned with hill-climbing (`bnlearn::hc`) scored by Gaussian BIC:

| DAG | Hypothesis | Parents of stress |
|---|---|---|
| **Learned** | Standard tier ordering, no extra rules | 9 variables |
| **Mediated** | School conditions affect stress *only* through psychological state | Psychological variables only |
| **Reverse** | Stress comes *before* anxiety, depression, sleep, and self-esteem | No psychological variables |
| **Mixed** | Sleep and self-esteem cause stress; anxiety and depression result from it | No anxiety or depression |

**Result:** the **Reverse DAG fit best**. The ranking was **Reverse > Mixed > Learned > Mediated**, and it held up under three checks:
- re-learning each DAG with 20 random restarts, to avoid getting stuck in a local optimum
- a custom **mixed-type BIC** that models each variable by its real type (logistic for binary, ordinal logistic for stress, linear otherwise)
- an unconstrained search, which also pointed stress → anxiety, depression, and sleep. It drew 37 arrows that break the tier ordering, including 12 that are impossible in time, which confirms BIC alone can't settle causal direction.

### 3. Identification
Using `dagitty` on the Reverse DAG:
- Every school factor has directed causal paths to stress.
- **Backdoor adjustment sets** were found for each factor's **total effect** (controls for confounders) and **direct effect** (also controls for mediators).

### 4. Estimation
- **Naive, total, direct, and indirect effects** from linear regression, where indirect = total − direct. These are reported both in raw units and standardized.
- **Ordinal logistic regression** (`MASS::polr`) as a robustness check, since stress is ordered categories rather than a continuous number. Effects are turned into probability changes with `marginaleffects`.

| Factor | Odds ratio | Change in P(high stress) |
|---|---|---|
| Study load | 1.77 | +3.8 pp |
| Career concerns | 1.64 | +3.3 pp |
| Teacher relationship | 0.42 | −10.1 pp |

### 5. Sensitivity Analysis
Three questions test whether the conclusions depend on choices made along the way:

| Question | Method | Result |
|---|---|---|
| Does it depend on how stress is measured? | Ordinal vs. linear models | Same signs and significance |
| Does it depend on the chosen DAG? | Re-find the adjustment sets and re-estimate under 4 alternative DAGs | Study load and career concerns are identical under every DAG. The teacher effect ranges from −0.09 to −0.17, because the data can't tell whether teacher relationships affect bullying or the reverse (ΔBIC = 0.16). |
| Does it depend on unmeasured variables? | `sensemakr`, benchmarked against hidden confounders 1–3× as strong as mental health history | Study load, career concerns, and the conservative teacher estimate survive a confounder 3× as strong. The larger teacher estimate (−0.17) loses significance at 2×. |

---

## Limitations

- **Cross-sectional self-report data:** the time-ordered tiers come from theory, not from measuring students at different points in time.
- Some arrows the unconstrained search drew (e.g., stress → academic performance, or psychological state → *how students perceive* school conditions) suggest reverse causation or perception bias.
- The direction between teacher relationships and bullying or peer pressure can't be determined from the data, so the teacher effect is reported as a range.
- The main estimates treat the ordinal outcome as continuous. The ordinal model confirms the direction and significance of the results.

---
