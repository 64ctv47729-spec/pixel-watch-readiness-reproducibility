#!/usr/bin/env python3
"""Outcome-specific inverse-probability-weighted Aim 1 sensitivity.

For each ordinal outcome, model availability of the transcribed Readiness
exposure from protocol time, episode, and that observed wellness outcome.
Stabilised weights use time/episode in the numerator and the full model in the
denominator, are truncated at the 1st/99th percentiles, and enter an explicitly
labelled weighted pseudo-likelihood cumulative-logit random-intercept model.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm

from independent_models import OUTCOME_LEVELS, bh, fit_clmm, prepare

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
MOD = OUT / "models"


def main() -> None:
    ds = pd.read_csv(OUT / "analytic_dataset.csv", parse_dates=["q_date"])
    ds.loc[ds.Stress == "Very Stressed", "Stress"] = "Stressed"
    ds.loc[ds.Body == "Very Sore/Tight", "Body"] = "Sore/Tight"
    ds["avail"] = ds.Google_Readiness.notna().astype(int)
    ds["day7"] = ds.protocol_day / 7.0
    ds["second_episode"] = (ds.episode > 1).astype(int)

    rows, weight_rows = [], []
    for outcome, levels in OUTCOME_LEVELS.items():
        work = ds.loc[ds[outcome].notna()].copy()
        work["outcome_score"] = pd.Categorical(
            work[outcome], categories=levels, ordered=True
        ).codes
        work = work.loc[work.outcome_score >= 0].copy()
        numerator = smf.gee(
            "avail ~ day7 + second_episode", groups="Participant_ID", data=work,
            family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()
        ).fit()
        denominator = smf.gee(
            "avail ~ day7 + second_episode + outcome_score", groups="Participant_ID",
            data=work, family=sm.families.Binomial(),
            cov_struct=sm.cov_struct.Exchangeable()
        ).fit()
        selected = work.loc[work.avail == 1].copy()
        pn = np.clip(numerator.predict(selected), 1e-4, 1 - 1e-4)
        pdn = np.clip(denominator.predict(selected), 1e-4, 1 - 1e-4)
        selected["ipw"] = pn / pdn
        qlo, qhi = selected.ipw.quantile([0.01, 0.99])
        selected["ipw"] = selected.ipw.clip(qlo, qhi)
        weight_rows.extend({
            "outcome": outcome, "Participant_ID": int(r.Participant_ID),
            "q_date": str(r.q_date.date()), "ipw": float(r.ipw),
        } for _, r in selected.iterrows())

        dat = prepare(selected)
        codes = pd.Categorical(dat[outcome], categories=levels, ordered=True).codes
        ok = codes >= 0
        dat = dat.loc[ok].copy()
        codes = codes[ok]
        x = dat[["r_within", "r_between", "day_std"]].to_numpy(float)
        fit = fit_clmm(codes, x, dat.Participant_ID.to_numpy(), nodes=30,
                       weights=dat.ipw.to_numpy())
        for j, term in enumerate(["r_within", "r_between", "day_std"]):
            beta, se = float(fit["beta"][j]), float(fit["se"][j])
            p = 2 * norm.sf(abs(beta / se))
            rows.append({
                "analysis": "selection sensitivity: outcome-specific stabilised IPW pseudo-likelihood",
                "outcome": outcome, "term": term, "n": len(dat),
                "participants": dat.Participant_ID.nunique(), "beta": beta,
                "se_naive_weighted": se, "OR": math.exp(beta),
                "lo_naive_weighted": math.exp(beta - 1.96 * se),
                "hi_naive_weighted": math.exp(beta + 1.96 * se), "p_naive_weighted": p,
                "weight_min": float(dat.ipw.min()), "weight_max": float(dat.ipw.max()),
                "converged": fit["success"], "max_grad": fit["max_grad"],
                "note": "point-estimate sensitivity; Hessian CI is not a robust IPW CI",
            })
    result = pd.DataFrame(rows)
    for term in ["r_within", "r_between"]:
        ii = result.term == term
        result.loc[ii, "q_naive_weighted"] = bh(result.loc[ii, "p_naive_weighted"].to_numpy())
    result.to_csv(MOD / "aim1_ipw_selection_sensitivity.csv", index=False)
    pd.DataFrame(weight_rows).to_csv(MOD / "aim1_ipw_selection_weights.csv", index=False)
    print(result.loc[result.term == "r_within"].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
