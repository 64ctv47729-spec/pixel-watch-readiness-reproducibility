#!/usr/bin/env python3
"""Additional analyses requested during peer review.

These outputs distinguish questionnaire-day transcription from device-score
availability, describe the calibration-compatible observation windows, and
compare wellness responses on days with and without a transcribed Readiness
score.  Blank classifications are operational descriptions, not inferred
causes of missingness.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.genmod.cov_struct import Independence
from statsmodels.genmod.generalized_estimating_equations import OrdinalGEE


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
REVIEW = OUT / "reviewer"
REVIEW.mkdir(parents=True, exist_ok=True)

LEVELS = {
    "Fatigue": ["Very Tired", "Tired", "Normal", "Fresh", "Very Fresh"],
    "Stress": ["Stressed", "Normal", "Relaxed", "Very Relaxed"],
    "Body": ["Sore/Tight", "Normal", "Feeling Good", "Feeling Great"],
    "Sleep": ["Restless Sleep", "Sleep Onset/Early Waking", "Restful", "Very Restful"],
}


def n_pct(frame: pd.DataFrame, label: str) -> dict[str, object]:
    n = int(len(frame))
    scored = int(frame["readiness_available"].sum())
    return {
        "window": label,
        "questionnaire_days": n,
        "readiness_recorded": scored,
        "readiness_blank": n - scored,
        "recorded_pct": 100 * scored / n if n else np.nan,
    }


def main() -> None:
    ds = pd.read_csv(OUT / "analytic_dataset.csv", parse_dates=["q_date"])
    ds["Stress"] = ds["Stress"].replace({"Very Stressed": "Stressed"})
    ds["Body"] = ds["Body"].replace({"Very Sore/Tight": "Sore/Tight"})
    ds["readiness_available"] = ds["Google_Readiness"].notna().astype(int)
    ds["episode_key"] = ds["Participant_ID"].astype(str) + "-" + ds["episode"].astype(str)

    # Denominators that answer different availability questions.
    rows = [
        n_pct(ds, "All questionnaire-days"),
        n_pct(ds[ds.protocol_day.between(1, 7)], "Protocol days 1-7"),
        n_pct(ds[ds.protocol_day.between(8, 14)], "Protocol days 8-14"),
        n_pct(ds[ds.protocol_day >= 8], "Protocol day 8 onward"),
        n_pct(ds[(ds.episode == 1) & ds.protocol_day.between(8, 14)],
              "Episode 1, protocol days 8-14"),
    ]
    first_day = (
        ds.loc[ds.readiness_available.eq(1)]
        .groupby("episode_key")["protocol_day"].min()
        .rename("first_recorded_day")
    )
    ds = ds.join(first_day, on="episode_key")
    after_first = ds[
        ds.first_recorded_day.notna() & (ds.protocol_day >= ds.first_recorded_day)
    ]
    rows.append(n_pct(after_first, "At or after first recorded score within score-yielding episode"))
    pd.DataFrame(rows).to_csv(REVIEW / "readiness_availability_windows.csv", index=False)

    # Classify blanks relative to the first observed score.  This deliberately
    # avoids assigning an undocumented technical or behavioural cause.
    def blank_class(row: pd.Series) -> str:
        if row.readiness_available:
            return "Readiness recorded"
        if pd.isna(row.first_recorded_day):
            return "Blank in episode with no recorded score"
        if row.protocol_day < row.first_recorded_day:
            return "Blank before first recorded score"
        return "Blank after first recorded score"

    ds["availability_class"] = ds.apply(blank_class, axis=1)
    class_order = [
        "Readiness recorded",
        "Blank before first recorded score",
        "Blank in episode with no recorded score",
        "Blank after first recorded score",
    ]
    cls = (
        ds.groupby("availability_class", observed=True)
        .agg(questionnaire_days=("availability_class", "size"),
             participants=("Participant_ID", "nunique"),
             episodes=("episode_key", "nunique"))
        .reindex(class_order).reset_index()
    )
    cls["pct_all_questionnaire_days"] = 100 * cls.questionnaire_days / len(ds)
    cls.to_csv(REVIEW / "readiness_blank_classification.csv", index=False)

    first = first_day.reset_index()
    first["first_day_group"] = pd.cut(
        first.first_recorded_day,
        bins=[0, 1, 7, 14, np.inf],
        labels=["Day 1", "Days 2-7", "Days 8-14", "Day 15 onward"],
    )
    first_dist = (
        first.groupby("first_day_group", observed=False).size()
        .rename("episodes").reset_index()
    )
    first_dist["pct_score_yielding_episodes"] = 100 * first_dist.episodes / len(first)
    first_dist.to_csv(REVIEW / "readiness_first_score_distribution.csv", index=False)

    per_person = (
        ds.groupby("Participant_ID")
        .agg(questionnaire_days=("readiness_available", "size"),
             readiness_days=("readiness_available", "sum"))
        .reset_index()
    )
    per_person["readiness_pct"] = 100 * per_person.readiness_days / per_person.questionnaire_days
    per_person.to_csv(REVIEW / "readiness_days_per_participant.csv", index=False)
    freq = (
        per_person.groupby("readiness_days").size().rename("participants").reset_index()
    )
    freq["pct_participants"] = 100 * freq.participants / len(per_person)
    freq.to_csv(REVIEW / "readiness_days_per_participant_distribution.csv", index=False)

    summary = {
        "participants": int(len(per_person)),
        "median_readiness_days": float(per_person.readiness_days.median()),
        "q1_readiness_days": float(per_person.readiness_days.quantile(.25)),
        "q3_readiness_days": float(per_person.readiness_days.quantile(.75)),
        "minimum_readiness_days": int(per_person.readiness_days.min()),
        "maximum_readiness_days": int(per_person.readiness_days.max()),
        "participants_with_zero_readiness_days": int(per_person.readiness_days.eq(0).sum()),
        "median_readiness_pct": float(per_person.readiness_pct.median()),
        "q1_readiness_pct": float(per_person.readiness_pct.quantile(.25)),
        "q3_readiness_pct": float(per_person.readiness_pct.quantile(.75)),
        "score_yielding_episodes": int(first_day.size),
        "episodes_without_recorded_score": int(ds.episode_key.nunique() - first_day.size),
    }
    (REVIEW / "readiness_availability_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    # Category counts are the most transparent unadjusted comparison.
    cat_rows: list[dict[str, object]] = []
    model_rows: list[dict[str, object]] = []
    for outcome, levels in LEVELS.items():
        d = ds[ds[outcome].isin(levels)].copy()
        for available, label in [(0, "Readiness blank"), (1, "Readiness recorded")]:
            z = d[d.readiness_available.eq(available)]
            denom = len(z)
            for level in levels:
                count = int(z[outcome].eq(level).sum())
                cat_rows.append({
                    "outcome": outcome,
                    "availability": label,
                    "category": level,
                    "n": count,
                    "denominator": denom,
                    "pct": 100 * count / denom if denom else np.nan,
                })

        d["y"] = pd.Categorical(d[outcome], categories=levels, ordered=True).codes
        d["protocol_week"] = d.protocol_day / 7.0
        d["second_episode"] = (d.episode > 1).astype(int)
        # OrdinalGEE creates threshold intercepts internally; a separate
        # constant is therefore intentionally omitted.
        model = OrdinalGEE.from_formula(
            "y ~ 0 + readiness_available + protocol_week + second_episode",
            groups="Participant_ID",
            data=d,
            cov_struct=Independence(),
        )
        result = model.fit(maxiter=200)
        term = "readiness_available"
        beta = float(result.params[term])
        se = float(result.bse[term])
        model_rows.append({
            "outcome": outcome,
            "n": int(len(d)),
            "participants": int(d.Participant_ID.nunique()),
            "term": term,
            "beta_log_odds_more_favourable": beta,
            "robust_se": se,
            "OR_more_favourable": float(np.exp(beta)),
            "ci_low": float(np.exp(beta - 1.96 * se)),
            "ci_high": float(np.exp(beta + 1.96 * se)),
            "p": float(result.pvalues[term]),
            "adjusted_for": "protocol day (per 7 days), second episode",
            "working_correlation": "independence; participant-cluster robust SE",
        })

    pd.DataFrame(cat_rows).to_csv(REVIEW / "wellness_by_readiness_availability.csv", index=False)
    pd.DataFrame(model_rows).to_csv(REVIEW / "wellness_availability_ordinal_gee.csv", index=False)

    print(json.dumps(summary, indent=2))
    print(pd.DataFrame(rows).round(2).to_string(index=False))
    print(pd.DataFrame(model_rows).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
