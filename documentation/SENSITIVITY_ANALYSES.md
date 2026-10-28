# Sensitivity analyses

## Aim 1 wellness

Within-person odds ratios per 10-point Readiness increment:

| Specification | Fatigue | Stress | Body | Sleep |
|---|---:|---:|---:|---:|
| Primary episode-centred | 1.192 (q=0.042) | 1.261 (q=0.014) | 1.114 (q=0.172) | 1.158 (q=0.085) |
| Episode 1, days 1–14 | 1.264 (q=0.037) | 1.248 (q=0.037) | 1.221 (q=0.046) | 1.099 (q=0.325) |
| Readiness=0 excluded | 1.231 (q=0.017) | 1.281 (q=0.010) | 1.093 (q=0.277) | 1.150 (q=0.124) |
| Participant-level centring | 1.195 (q=0.035) | 1.255 (q=0.014) | 1.127 (q=0.126) | 1.179 (q=0.047) |
| All duplicate-adjudicated dates excluded | 1.171 (q=0.079) | 1.234 (q=0.035) | 1.110 (q=0.197) | 1.156 (q=0.097) |
| Outcome-specific stabilised IPW | 1.193 (q=0.043) | 1.257 (q=0.016) | 1.117 (q=0.160) | 1.157 (q=0.086) |

The effect directions and magnitudes are broadly stable. Which individual wellness domains cross an FDR threshold is not: Body does so in the protocol restriction, Sleep under participant centring, and Fatigue does not after removing every duplicate-adjudicated date. The IPW point-estimate sensitivity closely matches the primary model.

## Aim 2 physiology

The protocol-restricted sample contains 170 participant-days from 30 identifiers. Independent estimates remain directionally consistent: RMSSD +7.44%, PNS +0.166, SNS −0.149 and heart rate −1.080 beats/min. The SD2 bootstrap interval includes zero in the independent rerun, even though its Wald/FDR test remains small; this reinforces interpreting the HRV measures as a correlated pattern.

Timing/linkage and quality sensitivities for RMSSD and heart rate:

| Restriction | n (identifiers) | RMSSD effect | Heart-rate effect |
|---|---:|---:|---:|
| Same local date, primary | 212 (31) | +6.20% | −0.99 beats/min |
| Artefact correction ≤5% | 207 (31) | +6.05% | −1.00 beats/min |
| Exclude duplicate-adjudicated dates | 199 (31) | +6.03% | −0.97 beats/min |
| Questionnaire after Kubios | 209 (31) | +6.24% | −1.04 beats/min |
| After and within 4 hours | 172 (26) | +7.55% | −1.19 beats/min |
| After and within 2 hours | 161 (25) | +8.40% | −1.20 beats/min |
| After and within 30 minutes | 123 (22) | +8.41% | −1.32 beats/min |

All displayed timing-window intervals exclude the null. Full results for all eight outcomes are in `output/models/aim2_timing_window_sensitivity.csv`; Readiness=0 and Kubios-quality sensitivities are in `output/models/aim2_independent_all.csv`.

## Aim 3 prediction

The submitted row-weighted point estimates exactly reproduce. After matching the cluster-bootstrap calculation to that same row-weighted estimand, Fatigue, Stress and Body remain non-significant and Sleep remains slightly worse. Participant-weighted sensitivity gives the same direction and conclusion for Sleep (difference +0.0021, 95% CI +0.0004 to +0.0042, FDR q=0.036).

The exploratory forward temporal analysis within known participants also shows no consistent benefit: RPS differences are +0.0044 (Fatigue), −0.0011 (Stress), +0.0117 (Body) and +0.0007 (Sleep).

Adding Google Sleep Score to the base model before Daily Readiness did not reveal incremental predictive value from Readiness. Across five-fold participant-grouped cross-validation repeated 20 times, the expanded-minus-base RPS differences were +0.0020 for Fatigue (95% CI +0.0004 to +0.0038; q=0.032), +0.0009 for Stress (−0.0035 to +0.0056; q=0.806), −0.0003 for Body (−0.0024 to +0.0018; q=0.806), and +0.0011 for Sleep (+0.0004 to +0.0019; q=0.002). Negative values favour Readiness. No validation fold failed. Watch-derived activity data were unavailable; self-reported exercise remains in the primary base model.

## Readiness availability and recorded-versus-blank days

Daily Readiness was transcribed on 220 of 410 questionnaire-days (53.7%). Availability was 123/263 (46.8%) on protocol days 1-7, 76/97 (78.4%) on protocol days 8-14, and 220/250 (88.0%) at or after the first recorded score within the 36 episodes that yielded a score. These are questionnaire-transcription proportions, not direct estimates of device-generation success. The 190 blank records comprise 96 before an episode's first score, 64 in eight episodes that never yielded a score, and 30 after an episode's first score; device-side causes cannot be separated.

Participants contributed a median of four recorded-score days (IQR 1-9; range 0-16); seven contributed none. In the full logistic GEE, the odds of a recorded score increased with protocol time (OR 1.62 per seven days, 95% CI 1.00-2.62; p=0.051), while the measured wellness variables were not clearly associated with recording. Direct adjusted ordinal-GEE contrasts likewise found no clear recorded-versus-blank difference for Fatigue (OR 1.26, 95% CI 0.82-1.94), Stress (1.06, 0.62-1.82), Body (1.56, 0.89-2.76) or Sleep (1.13, 0.75-1.69). These results reduce, but cannot eliminate, concern about informative missingness.

## Missingness limits

No tested sensitivity can identify the unobserved Readiness values or restore the seven identifiers with no score. These analyses demonstrate robustness to observed timing, protocol and selection variables, not proof that missingness is ignorable.
