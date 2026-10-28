#!/usr/bin/env python3
"""Aim 3 sensitivity: add Daily Readiness after Google Sleep Score.

This repeats the participant-grouped 5-fold x 20-repeat design on identical
complete rows.  Scaling is learned only from each training fold.  The base
model contains protocol day, physiology, behavioural covariates, and Google
Sleep Score; the expanded model additionally contains Daily Readiness.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from aim3_corrected import (
    BASE, BIN, BOOT, BOOT_SEED, CONT, LEVELS, N_FOLDS, SEEDS,
    fit_po, load, log_loss_rows, mae_rows, predict_proba, rps,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
CV = OUT / "cv_sleep_sensitivity"
CV.mkdir(parents=True, exist_ok=True)
MODEL_OUT = OUT / "models"


def main() -> None:
    ds = load()
    predictions: list[pd.DataFrame] = []
    folds: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []

    for outcome, levels in LEVELS.items():
        d = ds.dropna(subset=[outcome] + BASE + ["Google_sleep_score", "Google_Readiness"]).copy()
        d = d[d[outcome].isin(levels)].reset_index(drop=True)
        d["y"] = pd.Categorical(d[outcome], categories=levels, ordered=True).codes
        pids = np.array(sorted(d.Participant_ID.unique()))
        n_levels = len(levels)

        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            perm = rng.permutation(pids)
            fold_of = {pid: i % N_FOLDS for i, pid in enumerate(perm)}
            d["_fold"] = d.Participant_ID.map(fold_of)
            folds.extend({
                "outcome": outcome, "seed": seed,
                "Participant_ID": pid, "fold": fold_of[pid],
            } for pid in pids)

            for fold in range(N_FOLDS):
                train = d[d._fold.ne(fold)]
                test = d[d._fold.eq(fold)]
                if len(test) == 0 or train.y.nunique() < 2:
                    continue

                cont_cols = CONT + ["Google_sleep_score"]
                mu = train[cont_cols].mean()
                sd = train[cont_cols].std().replace(0, 1)
                readiness_mu = train.Google_Readiness.mean()
                readiness_sd = train.Google_Readiness.std()

                def design(frame: pd.DataFrame, expanded: bool) -> np.ndarray:
                    xc = (frame[cont_cols] - mu) / sd
                    x = pd.concat([xc, frame[BIN]], axis=1)
                    if expanded:
                        x["Google_Readiness_std"] = (
                            frame.Google_Readiness - readiness_mu
                        ) / readiness_sd
                    return x.to_numpy(float)

                seen = train.y.unique()
                base_fit = fit_po(train.y.to_numpy(), design(train, False), n_levels)
                expanded_fit = fit_po(train.y.to_numpy(), design(train, True), n_levels)
                if base_fit is None or expanded_fit is None:
                    failures.append({
                        "outcome": outcome, "seed": seed, "fold": fold,
                        "base_failed": base_fit is None,
                        "expanded_failed": expanded_fit is None,
                        "n_train": len(train), "n_test": len(test),
                    })
                    continue

                y = test.y.to_numpy()
                pb = predict_proba(base_fit, design(test, False), n_levels, seen)
                pe = predict_proba(expanded_fit, design(test, True), n_levels, seen)
                predictions.append(pd.DataFrame({
                    "outcome": outcome,
                    "seed": seed,
                    "fold": fold,
                    "Participant_ID": test.Participant_ID.to_numpy(),
                    "q_date": test.q_date.dt.date.astype(str).to_numpy(),
                    "y": y,
                    "rps_base_sleep": rps(pb, y),
                    "rps_exp_sleep_readiness": rps(pe, y),
                    "ll_base_sleep": log_loss_rows(pb, y),
                    "ll_exp_sleep_readiness": log_loss_rows(pe, y),
                    "mae_base_sleep": mae_rows(pb, y),
                    "mae_exp_sleep_readiness": mae_rows(pe, y),
                }))

    pred = pd.concat(predictions, ignore_index=True)
    pred.to_csv(CV / "row_predictions.csv.gz", index=False, compression="gzip")
    pd.DataFrame(folds).to_csv(CV / "fold_assignments.csv", index=False)
    pd.DataFrame(failures).to_csv(CV / "failed_folds.csv", index=False)

    result_rows: list[dict[str, object]] = []
    boot_rows: list[dict[str, object]] = []
    for index, outcome in enumerate(LEVELS):
        p = pred[pred.outcome.eq(outcome)]
        metric_cols = [
            "rps_base_sleep", "rps_exp_sleep_readiness",
            "ll_base_sleep", "ll_exp_sleep_readiness",
            "mae_base_sleep", "mae_exp_sleep_readiness",
        ]
        repeated = p.groupby("seed")[metric_cols].mean()
        row_mean = p.groupby(["Participant_ID", "q_date"])[
            ["rps_base_sleep", "rps_exp_sleep_readiness"]
        ].mean()
        row_mean["diff"] = row_mean.rps_exp_sleep_readiness - row_mean.rps_base_sleep
        clusters = row_mean.reset_index().groupby("Participant_ID").agg(
            diff_sum=("diff", "sum"),
            n_rows=("diff", "size"),
            base_mean=("rps_base_sleep", "mean"),
            exp_mean=("rps_exp_sleep_readiness", "mean"),
        )
        pids = clusters.index.to_numpy()
        rng = np.random.default_rng(BOOT_SEED + 100 + index)
        row_boot = np.empty(BOOT)
        participant_boot = np.empty(BOOT)
        for draw in range(BOOT):
            choice = rng.choice(pids, len(pids), replace=True)
            selected = clusters.loc[choice]
            row_boot[draw] = selected.diff_sum.sum() / selected.n_rows.sum()
            participant_boot[draw] = (selected.exp_mean - selected.base_mean).mean()
            boot_rows.append({
                "outcome": outcome, "draw": draw + 1,
                "row_weighted_diff": row_boot[draw],
                "participant_weighted_diff": participant_boot[draw],
            })

        base = repeated.rps_base_sleep.mean()
        expanded = repeated.rps_exp_sleep_readiness.mean()
        pw = row_mean.groupby("Participant_ID")[[
            "rps_base_sleep", "rps_exp_sleep_readiness"
        ]].mean().mean()
        p_boot = 2 * min((row_boot <= 0).mean(), (row_boot >= 0).mean())
        pw_p_boot = 2 * min(
            (participant_boot <= 0).mean(), (participant_boot >= 0).mean()
        )
        result_rows.append({
            "outcome": outcome,
            "n_rows": int(row_mean.shape[0]),
            "n_participants": int(p.Participant_ID.nunique()),
            "rps_base_sleep": base,
            "rps_exp_sleep_readiness": expanded,
            "rps_diff": expanded - base,
            "diff_boot_lo": float(np.quantile(row_boot, .025)),
            "diff_boot_hi": float(np.quantile(row_boot, .975)),
            "p_boot": max(float(p_boot), 1 / BOOT),
            "ll_base_sleep": repeated.ll_base_sleep.mean(),
            "ll_exp_sleep_readiness": repeated.ll_exp_sleep_readiness.mean(),
            "mae_base_sleep": repeated.mae_base_sleep.mean(),
            "mae_exp_sleep_readiness": repeated.mae_exp_sleep_readiness.mean(),
            "pw_rps_base_sleep": pw.rps_base_sleep,
            "pw_rps_exp_sleep_readiness": pw.rps_exp_sleep_readiness,
            "pw_rps_diff": pw.rps_exp_sleep_readiness - pw.rps_base_sleep,
            "pw_diff_boot_lo": float(np.quantile(participant_boot, .025)),
            "pw_diff_boot_hi": float(np.quantile(participant_boot, .975)),
            "pw_p_boot": max(float(pw_p_boot), 1 / BOOT),
        })

    results = pd.DataFrame(result_rows)
    results["q_boot"] = multipletests(results.p_boot, method="fdr_bh")[1]
    results["pw_q_boot"] = multipletests(results.pw_p_boot, method="fdr_bh")[1]
    results.to_csv(MODEL_OUT / "aim3_sleep_score_sensitivity.csv", index=False)
    pd.DataFrame(boot_rows).to_csv(
        CV / "cluster_bootstrap_draws.csv.gz", index=False, compression="gzip"
    )
    meta = {
        "n_folds": N_FOLDS,
        "n_repeats": len(SEEDS),
        "seeds": SEEDS,
        "bootstrap": BOOT,
        "bootstrap_seed_base": BOOT_SEED + 100,
        "failed_folds": len(failures),
        "base_model_addition": "standardized Google Sleep Score",
        "expanded_model_addition": "standardized Google Daily Readiness",
    }
    (CV / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    print(results.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
