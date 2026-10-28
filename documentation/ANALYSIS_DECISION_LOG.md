# Analysis decision log

## Status and timing

This decision log was fixed on 27 October 2026 after inventorying the available files and auditing the provenance of the submitted analysis, but before executing the new independent models in this audit directory. The analyst had access to the submitted manuscript and to earlier result files; therefore this is a locked reproducibility and appropriateness plan, not a claim of blinded preregistration. Model choices below follow the study questions, protocol, data structure and the audit brief rather than statistical significance.

## Source hierarchy and submitted target

The reporting target is the journal copy of `HRV_Readiness_Manuscript_Sensors_131026.docx` stored in `CERBERUS/Drafts/HRV Study 2026`, not the similarly named working copy in `final_analysis_202609`. Raw date-keyed questionnaire and Kubios exports take precedence over derived files. The consolidated workbook is used for Fitabase only where no full-cohort native export is available. Previous scripts and outputs are evidence of what was done, not evidence that the results are correct.

## Data construction

- Rebuild the participant-day dataset from the raw questionnaire export, raw Kubios export and the dated Fitabase block using participant identity and local calendar date in Europe/Dublin.
- Never join questionnaire and physiology by row position.
- Retain one questionnaire per participant-date using this deterministic order: most complete wellness record, Readiness present, earliest start time, lowest immutable respondent identifier.
- Define a new monitoring episode after a gap greater than 30 calendar days. Calculate protocol day as calendar day within episode, not as a sequential record index.
- Keep the value Readiness = 0 in the primary analysis and exclude it in sensitivity analysis because it may be a transcription sentinel.
- Do not impute Readiness.

## Primary and sensitivity samples

- Primary: all keyed participant-day records, using outcome-specific complete cases.
- Protocol-informed sensitivity: first episode only and protocol days 1 to 14. For Aim 2 this also requires an available questionnaire and Kubios record on the same local date. This sample is not described as fully linkage verified because Google Readiness was self-transcribed and has no independent timestamped export.
- Additional sensitivities: remove all participant-dates requiring duplicate adjudication; remove Readiness = 0; use participant-level rather than episode-level centring; require Kubios quality at least 90 where applicable.

## Exposure decomposition

Express Readiness per 10 points. Within each participant-episode, define the within-person component as the daily value minus the episode mean. Define the between-episode component as the episode mean minus the grand mean. Models include standardised calendar protocol day and a participant random intercept.

## Aim 1 wellness estimand and model

The primary estimand is the subject-specific within-person cumulative odds ratio for a more favourable ordinal wellness category per 10-point higher-than-usual Readiness. Fit cumulative-logit mixed models:

`ordered outcome ~ within Readiness + between Readiness + standardised protocol day + (1 | participant)`

Use the data-dictionary order. Merge an extreme category with its adjacent category only when the full rebuilt dataset contains fewer than five observations in that extreme category. Report odds ratios, Wald 95% confidence intervals, two-sided p-values, observations, participants, convergence diagnostics and random-intercept variation. Apply Benjamini-Hochberg correction across the four within-person tests; treat the four between-person tests as a separate family. Assess proportional odds with threshold-specific mixed logistic models. Sensitivities include ordinal marginal models, numeric-score mixed models and the samples above.

## Aim 2 physiology estimand and model

The primary estimand is the subject-specific within-person difference in independently derived physiology per 10-point higher-than-usual Readiness. Fit participant-random-intercept linear mixed models using restricted maximum likelihood. Natural-log transform RMSSD, SDNN, SD1 and SD2; report these as `100 * (exp(beta) - 1)` percent differences. Analyse PNS index, SNS index, heart rate and respiratory rate in native units. Use participant-cluster bootstrap percentile intervals with 2,000 resamples, alongside conventional mixed-model uncertainty. Apply Benjamini-Hochberg correction across the eight within-person outcomes and separately across the eight between-person outcomes. Treat RMSSD and SD1 as mathematically dependent signals. Kubios Readiness is supplementary context only.

## Aim 3 prediction target and validation

The main prediction target is same-day wellness in a previously unseen participant. Compare an ordinal base model containing protocol day, exercise, alcohol, illness, log RMSSD, heart rate and respiratory rate with the same model plus raw Readiness. Use identical complete cases and identical participant-grouped folds for both arms. Use five folds and 20 fixed-seed repeats. Estimate all scaling parameters from the training fold only. Save fold assignments, row-level class probabilities and failures. The primary metric is row-weighted ranked probability score; also report participant-weighted RPS, log loss and mean absolute ordinal error. Define difference as expanded minus base, so negative values favour Readiness. Use 2,000 participant-cluster bootstrap resamples and Benjamini-Hochberg correction across four outcomes. Treat forward prediction within known participants as exploratory.

## Missingness

Model the availability of a transcribed Readiness value with logistic GEE, exchangeable working correlation and participant-clustered robust standard errors. Use calendar protocol day, episode number and observed wellness predictors. Describe statistically detectable availability patterns as inconsistent with MCAR or potentially informative, never as proof of MNAR. A stabilised, truncated inverse-probability weighted Aim 1 analysis is a sensitivity only.

## Reproduction versus correction

First reproduce the scripts and values that support the submitted manuscript. Then independently refit the same primary specifications and protocol-informed sensitivities. Any discrepancy is classified by magnitude and scientific consequence; models are not changed to obtain statistical significance.

## Post-audit author decisions

These decisions were received after the audit models had been completed and therefore do not alter the locked analytical plan above.

- Treat the separately calculated demographic values in the manuscript as correct on author confirmation.
- Accept the author’s confirmation of 38 distinct eligible, consented participants, implemented recruitment criteria, administered wellness wording, questionnaire order, firmware and common-software claims.
- Use the 19 August 2025 protocol and manuscript where the earlier checklist conflicts.
- Retain all same-local-date records as the primary analytical target.
- Describe the episode-mean coefficient as between-episode rather than purely between-person.
- Use row-weighted RPS as the Aim 3 primary estimand and report participant-weighted results separately.
- Omit the later protocol’s greater-than-5% corrected-beat rule from the manuscript and primary inclusion criteria.
- Treat questionnaire delays as timestamp artefacts, retain “immediately after,” and remove the timestamp-derived delay statistics from the revision.
- Preserve associational and convergent-validity language.

