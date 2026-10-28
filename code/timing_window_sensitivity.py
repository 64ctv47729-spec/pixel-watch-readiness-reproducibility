#!/usr/bin/env python3
"""Aim 2 sensitivity to questionnaire--Kubios timing restrictions."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from independent_models import (LOG_OUTCOMES, RAW_OUTCOMES, bh,
                                fit_random_intercept_reml, prepare)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"


def main() -> None:
    ds = pd.read_csv(OUT / "analytic_dataset.csv", parse_dates=["q_date"])
    dup = pd.read_csv(OUT / "audit" / "duplicate_questionnaire_adjudication.csv")
    dup_keys = set(zip(dup.Participant_ID.astype(int), pd.to_datetime(dup.q_date).dt.date))
    not_adjudicated = pd.Series([
        (int(r.Participant_ID), r.q_date.date()) not in dup_keys for _, r in ds.iterrows()
    ], index=ds.index)
    windows = {
        "same local date (primary)": ds.RMSSD.notna(),
        "protocol QC: artefact <=5%": ds.RMSSD.notna() & ds.Artefact_pct.le(5),
        "exclude duplicate-adjudicated dates": ds.RMSSD.notna() & not_adjudicated,
        "questionnaire after Kubios": ds.RMSSD.notna() & (ds.q_minus_k_minutes >= 0),
        "after and within 4 hours": ds.RMSSD.notna() & ds.q_minus_k_minutes.between(0, 240),
        "after and within 2 hours": ds.RMSSD.notna() & ds.q_minus_k_minutes.between(0, 120),
        "after and within 30 minutes": ds.RMSSD.notna() & ds.q_minus_k_minutes.between(0, 30),
    }
    rows = []
    for label, keep in windows.items():
        dat = ds.loc[keep & ds.Google_Readiness.notna()].copy()
        dat = prepare(dat)
        for oc in LOG_OUTCOMES + RAW_OUTCOMES:
            yvar = oc
            if oc in LOG_OUTCOMES:
                dat[f"log_{oc}"] = np.log(dat[oc])
                yvar = f"log_{oc}"
            dd = dat.loc[dat[yvar].notna() & np.isfinite(dat[yvar])].copy()
            x = np.column_stack([np.ones(len(dd)),
                                 dd[["r_within", "r_between", "day_std"]].to_numpy(float)])
            fit = fit_random_intercept_reml(dd[yvar].to_numpy(float), x,
                                             dd.Participant_ID.to_numpy())
            beta, se = float(fit["beta"][1]), float(fit["se"][1])
            p = 2 * norm.sf(abs(beta / se))
            rows.append({
                "timing_window": label, "outcome": oc, "n": len(dd),
                "participants": dd.Participant_ID.nunique(), "beta": beta,
                "se_wald": se, "lo_wald": beta - 1.96 * se,
                "hi_wald": beta + 1.96 * se, "p_wald": p,
                "effect": 100 * (math.exp(beta) - 1) if oc in LOG_OUTCOMES else beta,
                "effect_lo_wald": 100 * (math.exp(beta - 1.96 * se) - 1)
                if oc in LOG_OUTCOMES else beta - 1.96 * se,
                "effect_hi_wald": 100 * (math.exp(beta + 1.96 * se) - 1)
                if oc in LOG_OUTCOMES else beta + 1.96 * se,
                "scale": "percent" if oc in LOG_OUTCOMES else "native units",
                "converged": fit["success"],
            })
    res = pd.DataFrame(rows)
    for label in res.timing_window.unique():
        ii = res.timing_window == label
        res.loc[ii, "q_wald"] = bh(res.loc[ii, "p_wald"].to_numpy())
    res.to_csv(OUT / "models" / "aim2_timing_window_sensitivity.csv", index=False)
    print(res.loc[res.outcome.isin(["RMSSD", "HR"])].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
