#!/usr/bin/env python3
"""Corrected Aim 3 reanalysis.

Does same-day Google Readiness improve OUT-OF-SAMPLE prediction of
ordinal wellness beyond behavioural + physiological covariates?

Main target: previously unseen participants (participant-grouped CV,
5 folds x 20 repeats, seeds 1..20, identical folds for the paired models).
Exploratory target: future days within known participants (temporally blocked).

Leakage controls: standardisation uses training-fold means/SDs only; category
merges were fixed a priori from the full-data audit (ANALYSIS_DECISION_LOG §4);
base and expanded models run on exactly the same rows; failed folds dropped
from BOTH arms and logged. Fold assignments and full row-level probability
predictions are saved. The primary cluster bootstrap uses the same row-weighted
estimand as the displayed base, expanded, and paired-difference estimates.
Participant-weighted estimates and bootstrap intervals are reported separately.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.miscmodels.ordinal_model import OrderedModel

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
CV = OUT / "cv_corrected"
CV.mkdir(parents=True, exist_ok=True)
CORRECTED = ROOT / "results" / "corrected" / "aim3"
CORRECTED.mkdir(parents=True, exist_ok=True)

N_FOLDS = 5
N_REPEATS = 20
SEEDS = list(range(1, N_REPEATS + 1))
BOOT = 2000
BOOT_SEED = 20261027

LEVELS = {
    "Fatigue": ["Very Tired", "Tired", "Normal", "Fresh", "Very Fresh"],
    "Stress": ["Stressed", "Normal", "Relaxed", "Very Relaxed"],
    "Body": ["Sore/Tight", "Normal", "Feeling Good", "Feeling Great"],
    "Sleep": ["Restless Sleep", "Sleep Onset/Early Waking", "Restful", "Very Restful"],
}
CONT = ["protocol_day", "logRMSSD", "HR", "Resp_Rate_breaths_per_min"]
BIN = ["Exercise_b", "Alcohol_b", "Illness_b"]
BASE = CONT + BIN
EXP = BASE + ["Google_Readiness_std"]


def load() -> pd.DataFrame:
    ds = pd.read_csv(OUT / "analytic_dataset.csv", parse_dates=["q_date"])
    ds["Stress"] = ds["Stress"].replace({"Very Stressed": "Stressed"})
    ds["Body"] = ds["Body"].replace({"Very Sore/Tight": "Sore/Tight"})
    ds["logRMSSD"] = np.log(ds["RMSSD"])
    for b, c in [("Exercise_b", "Exercise"), ("Alcohol_b", "Alcohol"), ("Illness_b", "Illness")]:
        ds[b] = ds[c].map({"Yes": 1, "No": 0})
    return ds


def fit_po(y: np.ndarray, X: np.ndarray, n_levels: int):
    try:
        m = OrderedModel(y, X, distr="logit")
        r = m.fit(method="bfgs", maxiter=400, disp=False)
        if not r.mle_retvals.get("converged", False):
            r2 = m.fit(method="nm", maxiter=4000, disp=False)
            if not r2.mle_retvals.get("converged", False):
                return None
            r = r2
        return (m, r)
    except Exception:
        return None


def predict_proba(fit, X: np.ndarray, n_levels: int, levels_seen: np.ndarray) -> np.ndarray:
    m, r = fit
    p_seen = r.predict(exog=X)          # columns = sorted(levels_seen)
    p = np.zeros((X.shape[0], n_levels))
    for j, lev in enumerate(np.sort(levels_seen)):
        p[:, int(lev)] = p_seen[:, j]
    return np.clip(p, 1e-12, 1)


def rps(p: np.ndarray, y: np.ndarray) -> np.ndarray:
    K = p.shape[1]
    cp = np.cumsum(p, axis=1)
    obs = (y[:, None] <= np.arange(K)[None, :]).astype(float)
    return ((cp - obs) ** 2)[:, :-1].sum(axis=1) / (K - 1)


def log_loss_rows(p: np.ndarray, y: np.ndarray) -> np.ndarray:
    return -np.log(p[np.arange(len(y)), y])


def mae_rows(p: np.ndarray, y: np.ndarray) -> np.ndarray:
    exp_cat = (p * np.arange(p.shape[1])[None, :]).sum(axis=1)
    return np.abs(exp_cat - y)


def run_grouped_cv(ds: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list, pd.DataFrame]:
    all_preds, fail_log, folds, boot_draws = [], [], [], []
    for oc, levels in LEVELS.items():
        d = ds.dropna(subset=[oc] + BASE + ["Google_Readiness"]).copy()
        d = d[d[oc].isin(levels)]
        d["y"] = pd.Categorical(d[oc], categories=levels, ordered=True).codes
        d = d.reset_index(drop=True)
        pids = np.array(sorted(d["Participant_ID"].unique()))
        K = len(levels)

        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            perm = rng.permutation(pids)
            fold_of = {pid: i % N_FOLDS for i, pid in enumerate(perm)}
            d["_fold"] = d["Participant_ID"].map(fold_of)
            folds.extend({"outcome": oc, "seed": seed,
                          "Participant_ID": pid, "fold": fold_of[pid]}
                         for pid in pids)
            for fold in range(N_FOLDS):
                tr = d[d["_fold"] != fold]
                te = d[d["_fold"] == fold]
                if len(te) == 0 or tr["y"].nunique() < 2:
                    continue
                mu, sd = tr[CONT].mean(), tr[CONT].std().replace(0, 1)
                r_mu, r_sd = tr["Google_Readiness"].mean(), tr["Google_Readiness"].std()

                def design(frame, expanded):
                    Xc = (frame[CONT] - mu) / sd
                    X = pd.concat([Xc, frame[BIN]], axis=1)
                    if expanded:
                        X["Google_Readiness_std"] = (frame["Google_Readiness"] - r_mu) / r_sd
                    return X.to_numpy(float)

                seen = tr["y"].unique()
                fit_b = fit_po(tr["y"].to_numpy(), design(tr, False), K)
                fit_e = fit_po(tr["y"].to_numpy(), design(tr, True), K)
                if fit_b is None or fit_e is None:
                    fail_log.append({"outcome": oc, "seed": seed, "fold": fold,
                                     "base_failed": fit_b is None, "exp_failed": fit_e is None,
                                     "n_train": len(tr), "n_test": len(te)})
                    continue
                pb = predict_proba(fit_b, design(te, False), K, seen)
                pe = predict_proba(fit_e, design(te, True), K, seen)
                y = te["y"].to_numpy()
                pred_frame = pd.DataFrame({
                    "outcome": oc, "seed": seed, "fold": fold,
                    "Participant_ID": te["Participant_ID"].to_numpy(),
                    "q_date": te["q_date"].dt.date.astype(str).to_numpy(),
                    "y": y,
                    "rps_base": rps(pb, y), "rps_exp": rps(pe, y),
                    "ll_base": log_loss_rows(pb, y), "ll_exp": log_loss_rows(pe, y),
                    "mae_base": mae_rows(pb, y), "mae_exp": mae_rows(pe, y),
                })
                for j in range(K):
                    pred_frame[f"base_p{j}"] = pb[:, j]
                    pred_frame[f"expanded_p{j}"] = pe[:, j]
                for j in range(K, 5):
                    pred_frame[f"base_p{j}"] = np.nan
                    pred_frame[f"expanded_p{j}"] = np.nan
                all_preds.append(pred_frame)

    preds = pd.concat(all_preds, ignore_index=True)
    preds.to_csv(CV / "grouped_cv_row_predictions.csv.gz", index=False, compression="gzip")
    fold_df = pd.DataFrame(folds)
    fold_df.to_csv(CV / "fold_assignments.csv", index=False)

    rows = []
    for oc in LEVELS:
        p = preds[preds["outcome"] == oc]
        # row-weighted: mean over rows within a repeat, then mean over repeats
        rep = p.groupby("seed")[["rps_base", "rps_exp", "ll_base", "ll_exp",
                                 "mae_base", "mae_exp"]].mean()
        # Average repeated predictions for each unique participant-date before
        # bootstrap so that repeats are Monte Carlo repetitions, not replicated
        # observations.
        row_mean = p.groupby(["Participant_ID", "q_date"])[["rps_base", "rps_exp"]].mean()
        row_mean["diff"] = row_mean["rps_exp"] - row_mean["rps_base"]
        cluster = row_mean.reset_index().groupby("Participant_ID").agg(
            diff_sum=("diff", "sum"), n_rows=("diff", "size"),
            base_mean=("rps_base", "mean"), exp_mean=("rps_exp", "mean"),
        )
        pids = cluster.index.to_numpy()
        rng = np.random.default_rng(BOOT_SEED + list(LEVELS).index(oc))
        row_boots = np.empty(BOOT)
        participant_boots = np.empty(BOOT)
        for b in range(BOOT):
            choice = rng.choice(pids, len(pids), replace=True)
            chosen = cluster.loc[choice]
            row_boots[b] = chosen["diff_sum"].sum() / chosen["n_rows"].sum()
            participant_boots[b] = (chosen["exp_mean"] - chosen["base_mean"]).mean()
            boot_draws.append({"outcome": oc, "draw": b + 1,
                               "row_weighted_diff": row_boots[b],
                               "participant_weighted_diff": participant_boots[b]})
        diff = rep["rps_exp"].mean() - rep["rps_base"].mean()
        p_boot = 2 * min((row_boots <= 0).mean(), (row_boots >= 0).mean())
        # participant-weighted (sensitivity): same convention for all three
        pw = row_mean.groupby("Participant_ID")[["rps_base", "rps_exp"]].mean().mean()
        p_pw = 2 * min((participant_boots <= 0).mean(),
                       (participant_boots >= 0).mean())
        rows.append({
            "outcome": oc,
            "n_rows": int(p.groupby(["Participant_ID", "q_date"]).ngroups),
            "n_participants": int(p["Participant_ID"].nunique()),
            "rps_base": rep["rps_base"].mean(), "rps_exp": rep["rps_exp"].mean(),
            "rps_diff": diff,
            "diff_boot_lo": np.quantile(row_boots, .025),
            "diff_boot_hi": np.quantile(row_boots, .975),
            "p_boot": max(p_boot, 1 / BOOT),
            "ll_base": rep["ll_base"].mean(), "ll_exp": rep["ll_exp"].mean(),
            "mae_base": rep["mae_base"].mean(), "mae_exp": rep["mae_exp"].mean(),
            "pw_rps_base": pw["rps_base"], "pw_rps_exp": pw["rps_exp"],
            "pw_rps_diff": pw["rps_exp"] - pw["rps_base"],
            "pw_diff_boot_lo": np.quantile(participant_boots, .025),
            "pw_diff_boot_hi": np.quantile(participant_boots, .975),
            "pw_p_boot": max(p_pw, 1 / BOOT),
            "weighting": "primary row-weighted; pw_* participant-weighted sensitivity",
        })
    res = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    res["q_boot"] = multipletests(res["p_boot"], method="fdr_bh")[1]
    res["pw_q_boot"] = multipletests(res["pw_p_boot"], method="fdr_bh")[1]
    return res, preds, fold_df, fail_log, pd.DataFrame(boot_draws)


def run_temporal_cv(ds: pd.DataFrame) -> pd.DataFrame:
    """Exploratory: future-day prediction within known participants.
    Forward chaining on protocol day within episode 1: train on days <= c,
    test on days > c, for cuts c in {7, 9, 11}. Same paired design."""
    rows = []
    for oc, levels in LEVELS.items():
        d = ds[(ds["episode"] == 1)].dropna(subset=[oc] + BASE + ["Google_Readiness"]).copy()
        d = d[d[oc].isin(levels)]
        d["y"] = pd.Categorical(d[oc], categories=levels, ordered=True).codes
        K = len(levels)
        per_cut = []
        for cut in (7, 9, 11):
            tr = d[d["protocol_day"] <= cut]
            te = d[(d["protocol_day"] > cut) &
                   (d["Participant_ID"].isin(tr["Participant_ID"]))]
            if len(te) < 15 or tr["y"].nunique() < 2:
                continue
            mu, sd = tr[CONT].mean(), tr[CONT].std().replace(0, 1)
            r_mu, r_sd = tr["Google_Readiness"].mean(), tr["Google_Readiness"].std()

            def design(frame, expanded):
                Xc = (frame[CONT] - mu) / sd
                X = pd.concat([Xc, frame[BIN]], axis=1)
                if expanded:
                    X["Google_Readiness_std"] = (frame["Google_Readiness"] - r_mu) / r_sd
                return X.to_numpy(float)

            seen = tr["y"].unique()
            fb = fit_po(tr["y"].to_numpy(), design(tr, False), K)
            fe = fit_po(tr["y"].to_numpy(), design(tr, True), K)
            if fb is None or fe is None:
                continue
            y = te["y"].to_numpy()
            per_cut.append({
                "cut": cut, "n_test": len(te),
                "rps_base": rps(predict_proba(fb, design(te, False), K, seen), y).mean(),
                "rps_exp": rps(predict_proba(fe, design(te, True), K, seen), y).mean(),
            })
        if per_cut:
            pc = pd.DataFrame(per_cut)
            rows.append({"outcome": oc, "cuts_used": len(pc),
                         "n_test_total": int(pc["n_test"].sum()),
                         "rps_base": float(np.average(pc["rps_base"], weights=pc["n_test"])),
                         "rps_exp": float(np.average(pc["rps_exp"], weights=pc["n_test"])),
                         "rps_diff": float(np.average(pc["rps_exp"] - pc["rps_base"],
                                                      weights=pc["n_test"]))})
    return pd.DataFrame(rows)


def main() -> None:
    ds = load()
    res, preds, fold_df, fail_log, boot_draws = run_grouped_cv(ds)
    res.to_csv(CORRECTED / "aim3_grouped_cv_corrected.csv", index=False)
    boot_draws.to_csv(CV / "cluster_bootstrap_draws.csv.gz", index=False,
                      compression="gzip")
    pd.DataFrame(fail_log).to_csv(CV / "failed_folds.csv", index=False)
    temporal = run_temporal_cv(ds)
    temporal.to_csv(CORRECTED / "aim3_temporal_exploratory.csv", index=False)
    meta = {
        "n_folds": N_FOLDS, "n_repeats": N_REPEATS, "seeds": SEEDS,
        "bootstrap": BOOT, "bootstrap_seed_base": BOOT_SEED,
        "primary_weighting": "row weighted",
        "sensitivity_weighting": "participant weighted",
        "failed_folds": len(fail_log),
        "total_fold_fits": int(preds.groupby(["outcome", "seed", "fold"]).ngroups),
    }
    (CV / "aim3_meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))
    print(res.round(4).to_string(index=False))
    print("\nTEMPORAL (exploratory):")
    print(temporal.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
