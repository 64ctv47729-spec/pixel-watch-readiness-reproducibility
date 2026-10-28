#!/usr/bin/env python3
"""Google Readiness availability analysis.

Logistic GEE (exchangeable working correlation, robust SEs, clustered on
participant identifier) of availability on protocol day (per 7 days), episode
number, and same-day wellness scores. Language constraint: MCAR implausible /
systematic availability; never 'MNAR proven'.
Also: prespecified IPW sensitivity for the Aim 1 CLMM (weights exported for R).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"

ds = pd.read_csv(OUT / "analytic_dataset.csv", parse_dates=["q_date"])
ds["avail"] = ds["Google_Readiness"].notna().astype(int)
ds["day7"] = ds["protocol_day"] / 7.0
ds["second_episode"] = (ds["episode"] > 1).astype(int)
for oc in ["Fatigue", "Stress", "Body", "Sleep"]:
    ds[oc.lower() + "_s"] = ds[oc + "_score"]

d = ds.dropna(subset=["fatigue_s", "stress_s", "body_s", "sleep_s"]).copy()

model = smf.gee("avail ~ day7 + second_episode + fatigue_s + stress_s + body_s + sleep_s",
                groups="Participant_ID", data=d, family=sm.families.Binomial(),
                cov_struct=sm.cov_struct.Exchangeable())
res = model.fit()
tab = pd.DataFrame({
    "term": res.params.index, "beta": res.params.values, "se": res.bse.values,
    "OR": np.exp(res.params.values),
    "lo": np.exp(res.params.values - 1.96 * res.bse.values),
    "hi": np.exp(res.params.values + 1.96 * res.bse.values),
    "p": res.pvalues.values,
})
tab.to_csv(OUT / "models" / "missingness_gee.csv", index=False)

# availability by protocol day band and episode
band = pd.cut(d["protocol_day"], [0, 3, 7, 10, 14, 100], labels=["1-3", "4-7", "8-10", "11-14", ">14"])
avail_band = d.groupby(band, observed=True)["avail"].agg(["mean", "size"]).reset_index()
avail_band.to_csv(OUT / "models" / "missingness_by_day_band.csv", index=False)
by_ep = d.groupby("episode")["avail"].agg(["mean", "size"]).reset_index()

# participant availability
pa = ds.groupby("Participant_ID")["avail"].agg(["sum", "size"])
pa["pct"] = (pa["sum"] / pa["size"] * 100).round(1)
pa.reset_index().to_csv(OUT / "models" / "missingness_by_participant.csv", index=False)

# ---- IPW weights for Aim 1 sensitivity (stabilised, truncated 1st/99th pct)
ps = smf.gee("avail ~ day7 + second_episode", groups="Participant_ID", data=ds,
             family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
p_avail = ps.predict(ds)
marg = ds["avail"].mean()
w = np.where(ds["avail"] == 1, marg / p_avail, np.nan)
w = pd.Series(w).clip(*np.nanquantile(w, [.01, .99]))
ds["ipw"] = w
ds.loc[ds["avail"] == 1, ["Participant_ID", "q_date", "ipw"]].to_csv(
    OUT / "models" / "aim1_ipw_weights.csv", index=False)

summary = {
    "gee_n": int(res.nobs), "gee_participants": int(d["Participant_ID"].nunique()),
    "availability_overall_pct": round(float(ds["avail"].mean() * 100), 1),
    "availability_by_episode": by_ep.to_dict("records"),
    "or_day7": [round(float(np.exp(res.params["day7"])), 3),
                round(float(np.exp(res.params["day7"] - 1.96 * res.bse["day7"])), 3),
                round(float(np.exp(res.params["day7"] + 1.96 * res.bse["day7"])), 3),
                round(float(res.pvalues["day7"]), 5)],
    "ipw_range": [round(float(np.nanmin(w)), 3), round(float(np.nanmax(w)), 3)],
}
(OUT / "models" / "missingness_summary.json").write_text(json.dumps(summary, indent=2))
print(tab.round(4).to_string(index=False))
print(json.dumps(summary, indent=2))
