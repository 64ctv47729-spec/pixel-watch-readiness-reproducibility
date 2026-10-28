# Statistical methods audit

## Overall assessment

The primary model families are appropriate for the three stated questions and the repeated-measures structure. Aim 1 and Aim 2 estimates were reproduced from the rebuilt raw-data dataset and independently verified with separate likelihood/REML implementations. Aim 3 used a sound participant-grouped validation design, but its submitted bootstrap interval did not use the same weighting convention as its displayed point estimate. Correcting that issue changes uncertainty values but not the conclusion.

## Aim 1: ordinal wellness associations

Model: cumulative-logit mixed model with participant random intercept, episode-centred within-person Readiness, episode-mean Readiness, and standardised protocol day.

Strengths:

- Preserves ordinal outcomes instead of treating them as continuous.
- Separates day-to-day deviations from typical Readiness.
- Accounts for repeated observations with a participant random intercept.
- Reports observation and participant counts, convergence diagnostics, random-intercept variation, and FDR-adjusted results.
- Independent 30-node Gauss–Hermite likelihood estimates agree with the submitted `ordinal::clmm` results to negligible numerical error.

Qualifications:

- The “between-person” component is actually an episode-mean component. With 44 episodes from 38 identifiers, it mixes differences between people with differences between episodes for the six repeat identifiers. “Between-episode/typical-score association” is the precise label.
- Proportional-odds diagnostics are incomplete as evidence: one extreme fatigue threshold is too sparse among Readiness-complete days, and threshold-specific slopes vary, especially for stress. The directions are mostly consistent but precision is limited. The submitted manuscript says the assumption was evaluated but does not report the diagnostic results, and no Supporting Information file was found.
- Domain-level significance is sensitivity-dependent. Under protocol days 1–14, Body becomes FDR-significant; under participant-level centring, Sleep becomes marginally FDR-significant; after excluding all duplicate-adjudicated dates, Fatigue no longer survives FDR. The overall pattern of modest positive wellness associations is stable, but the claim that only Fatigue and Stress are associated should be presented as primary-analysis-specific.
- The exact item labels in Methods do not match the raw questionnaire labels.

Primary independent estimates per 10-point higher-than-usual Readiness:

| Outcome | OR (95% CI) | p | FDR q |
|---|---:|---:|---:|
| Fatigue/freshness | 1.192 (1.027–1.383) | 0.021 | 0.042 |
| Stress/relaxation | 1.261 (1.079–1.472) | 0.003 | 0.014 |
| Body soreness/recovery | 1.114 (0.954–1.300) | 0.172 | 0.172 |
| Sleep quality | 1.158 (0.992–1.353) | 0.064 | 0.085 |

Outcome-specific stabilised IPW pseudo-likelihood sensitivities give ORs of 1.193, 1.257, 1.117 and 1.157, respectively. This supports the primary point estimates under observed time/outcome-dependent selection, with the stated limitation that robust IPW uncertainty and unobserved missingness are not resolved.

## Aim 2: physiological associations

Model: participant-random-intercept linear mixed model with the same within/episode-mean decomposition and protocol day. RMSSD, SDNN, SD1 and SD2 are analysed on the natural-log scale and back-transformed; other outcomes remain in native units. Participant-cluster bootstrap percentile intervals use 2,000 resamples. BH correction is applied across eight within-person outcomes and separately across eight between-episode outcomes.

Assessment:

- The model and transformation choices are appropriate.
- Direct profiled REML reproduces submitted point estimates and Wald uncertainty to numerical precision.
- An independent 2,000-participant bootstrap produces small Monte Carlo differences in interval endpoints but identical conclusions.
- RMSSD and SD1 are mathematically dependent, and all HRV/autonomic outcomes are correlated. BH control is valid under common positive dependence, but the results should be interpreted as a coherent autonomic pattern, not eight independent confirmations.
- Same-date matching is defensible but broad. Restricting to questionnaires after Kubios and within 4 h, 2 h, or 30 min strengthens rather than removes the RMSSD/heart-rate associations. Excluding duplicate-adjudicated dates and applying artefact ≤5% also retain the conclusions.
- The later protocol’s <4.5-minute duration rule cannot be checked because duration is not exported.

The independent primary estimates remain: RMSSD +6.20%, SDNN +3.75%, SD1 +6.21%, SD2 +3.05%, PNS index +0.147, SNS index −0.125, heart rate −0.989 beats/min, and respiratory rate −0.057 breaths/min. All except respiratory rate survive within-family FDR correction.

## Aim 3: incremental prediction

Submitted design: ordinal base model versus the same model plus raw Readiness; five participant-grouped folds; 20 repeats with fixed seeds 1–20; identical folds and complete cases for the paired arms; training-fold-only scaling; held-out participants; ranked probability score (RPS) as the primary metric. There were 400 paired fold fits and no failed folds.

These leakage controls are implemented correctly. The base and expanded models use identical rows, transformations are learned from training data, and participants never cross train/test folds. Category merges were fixed before modelling from the full-data audit; because they define the outcome scale rather than estimate a predictor, the practical leakage risk is small, but this should remain explicit.

Submitted inference defect:

- The displayed base RPS, expanded RPS and difference are row-weighted.
- The submitted bootstrap first averaged by participant and therefore generated a participant-weighted interval and p-value.
- Thus, the reported point estimate and its uncertainty targeted different estimands.

The corrected analysis resamples participant blocks but computes each bootstrap draw as total loss difference divided by total resampled rows. A separate participant-weighted sensitivity is reported. Full class probabilities, fold assignments, seeds, losses, failures and 8,000 bootstrap draws are saved.

Corrected row-weighted results:

| Outcome | Base RPS | +Readiness RPS | Difference (95% cluster-bootstrap CI) | FDR q |
|---|---:|---:|---:|---:|
| Fatigue | 0.1316 | 0.1340 | +0.0024 (−0.0002 to +0.0050) | 0.144 |
| Stress | 0.1333 | 0.1349 | +0.0016 (−0.0033 to +0.0073) | 0.596 |
| Body | 0.1355 | 0.1346 | −0.0009 (−0.0031 to +0.0012) | 0.491 |
| Sleep | 0.1601 | 0.1614 | +0.0014 (+0.0004 to +0.0024) | 0.020 |

No outcome shows a predictive gain. Sleep remains slightly worse after adding Readiness. The scientific conclusion is unchanged.

## Multiplicity, clustering and uncertainty

- BH families match the manuscript: four within and four episode-mean tests for Aim 1; eight within and eight episode-mean tests for Aim 2; four outcomes for Aim 3.
- Participant clustering is used for model random effects or GEE robust variance and for cluster bootstrap resampling.
- Aim 1 primary CIs are model-based Wald intervals; small-cluster limitations should be acknowledged.
- Aim 2 bootstrap intervals correctly resample and relabel participants.
- Aim 3 is now internally consistent after the row-weighted correction.

## Missingness and causal interpretation

All results are associational and conditional on a recorded, manually transcribed Readiness score. The missingness process is not proven MNAR, and no causal effect of Readiness can be inferred. The independent reanalysis supports the manuscript’s central physiological and predictive conclusions but not a claim of complete representativeness or causal validation.
