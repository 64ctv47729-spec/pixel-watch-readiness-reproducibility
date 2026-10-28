# Reproducibility verification report

## Status

**Passed — 28 October 2026.**

The package was executed from its disclosure-reduced analytic dataset using the included scripts. Automated validation confirmed agreement with the frozen manuscript-supporting results and found no direct identifiers, email-like strings, local absolute paths, real study dates or symlinks in the release candidate.

## Data invariants

- Participant-day rows: 410
- Pseudonymous participants: 38
- Monitoring episodes: 44
- Rows with a recorded Daily Readiness score: 220
- Participants with a recorded Daily Readiness score: 31
- Real calendar dates included: none
- Direct response identifiers included: none
- Raw participant crosswalk included: no

## Analyses executed

- Independent 30-node Gauss–Hermite cumulative-logit mixed models for Aim 1.
- Independent profiled-REML random-intercept models for Aim 2.
- Participant-cluster bootstrap intervals for Aim 2 with 2,000 draws.
- Score-availability GEE and recorded-versus-blank ordinal GEE analyses.
- Outcome-specific inverse-probability-weighted Aim 1 sensitivity analyses.
- Timing-window, duplicate-date, protocol-restricted and participant-centring sensitivities.
- Participant-grouped five-fold cross-validation with 20 repeats for Aim 3.
- Row-weighted and participant-weighted ranked probability scores.
- Participant-cluster bootstrap inference for Aim 3 with 2,000 draws per outcome.
- Forward-blocked exploratory prediction and Google Sleep Score sensitivity analyses.
- Regeneration of analysis tables and manuscript figures.

There were no failed grouped-cross-validation folds in the primary or Sleep Score sensitivity analyses.

## Numerical verification

The following regenerated files were compared column-by-column with frozen expected results using numeric tolerances of `rtol=1e-7` and `atol=1e-9`:

- Aim 1 primary model results;
- Aim 2 primary model results;
- corrected Aim 3 results;
- Aim 3 Google Sleep Score sensitivity;
- score-availability GEE;
- recorded-versus-blank wellness comparisons.

All comparisons passed. The participant-weighted fatigue difference in Table S7 regenerated as +0.0007572933, reported as +0.0008. The Sleep Score sensitivity stress difference in Table S9 regenerated as +0.0008533962, reported as +0.0009. Both are therefore confirmed as unrounded-estimate differences.

## Integrity controls

- Every released file is listed in `MANIFEST.sha256`.
- `validate_package.sh` verifies data invariants, privacy guards, key numerical outputs and every manifest checksum.
- The package includes no symlinks and does not depend on the private project directory at runtime after dependencies are installed.
- Raw-to-analytic construction and provenance remain preserved in the private dated audit; the public candidate contains only the downstream analytic workflow.

## Qualification

The original Aim 1 and Aim 2 primary fits used R 4.3.3. R was not available in the packaging environment, so those scripts were retained but not re-executed during this final packaging run. The independent Python likelihood/REML implementation was executed and had already been shown to agree with the R estimates to negligible numerical error.

This technical verification does not itself authorise public dissemination of participant-level data. Consent, ethics, institutional governance, author approval and licensing remain release-gating decisions.

