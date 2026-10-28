#!/usr/bin/env python3
"""Validate privacy guards, data invariants and key numerical outputs."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_COLUMNS = {
    "respondent_id", "q_start", "q_end", "Measured_time", "k_time_local", "f_time"
}
TEXT_SUFFIXES = {".csv", ".json", ".md", ".py", ".R", ".sh", ".txt", ".cff", ".yml", ".yaml"}


def compare_csv(actual: Path, expected: Path, key_columns: list[str]) -> None:
    observed = pd.read_csv(actual).sort_values(key_columns).reset_index(drop=True)
    target = pd.read_csv(expected).sort_values(key_columns).reset_index(drop=True)
    if list(observed.columns) != list(target.columns):
        raise AssertionError(f"Column mismatch: {actual.name}")
    if observed.shape != target.shape:
        raise AssertionError(f"Shape mismatch: {actual.name}: {observed.shape} != {target.shape}")
    for column in observed.columns:
        if pd.api.types.is_numeric_dtype(target[column]):
            if not np.allclose(observed[column], target[column], rtol=1e-7, atol=1e-9, equal_nan=True):
                delta = np.nanmax(np.abs(observed[column].to_numpy(float) - target[column].to_numpy(float)))
                raise AssertionError(f"Numeric mismatch: {actual.name}:{column}; max delta={delta}")
        else:
            left = observed[column].fillna("<NA>").astype(str)
            right = target[column].fillna("<NA>").astype(str)
            if not left.equals(right):
                raise AssertionError(f"Text mismatch: {actual.name}:{column}")


def validate_data() -> None:
    data = pd.read_csv(ROOT / "data/analytic_dataset.csv", parse_dates=["q_date"])
    assert len(data) == 410
    assert data.Participant_ID.nunique() == 38
    assert data.groupby("Participant_ID").episode.max().sum() == 44
    assert data.Google_Readiness.notna().sum() == 220
    assert not FORBIDDEN_COLUMNS.intersection(data.columns)
    assert data.Participant_ID.min() == 1 and data.Participant_ID.max() == 38
    assert not data.q_date.dt.year.isin([2024, 2025, 2026]).any()
    assert data.q_date.notna().all()


def validate_privacy_scan() -> None:
    email_pattern = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
    findings = []
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if any(part in {".git", ".venv", "__pycache__"} for part in relative.parts):
            continue
        if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        text = path.read_text(errors="ignore")
        if email_pattern.search(text):
            findings.append(f"email-like string in {relative}")
        if "/Users/" in text or "CloudStorage" in text:
            findings.append(f"local absolute path in {relative}")
    if findings:
        raise AssertionError("Privacy scan failed:\n" + "\n".join(findings))


def validate_outputs() -> None:
    required = [
        ROOT / "output/models/aim1_independent_primary.csv",
        ROOT / "output/models/aim2_independent_primary.csv",
        ROOT / "results/corrected/aim3/aim3_grouped_cv_corrected.csv",
        ROOT / "output/models/aim3_sleep_score_sensitivity.csv",
        ROOT / "output/models/missingness_gee.csv",
        ROOT / "output/reviewer/wellness_availability_ordinal_gee.csv",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise AssertionError("Missing generated outputs: " + ", ".join(missing))
    compare_csv(required[0], ROOT / "expected_results/aim1_primary.csv", ["outcome", "term"])
    compare_csv(required[1], ROOT / "expected_results/aim2_primary.csv", ["outcome", "term"])
    compare_csv(required[2], ROOT / "expected_results/aim3_primary.csv", ["outcome"])
    compare_csv(required[3], ROOT / "expected_results/aim3_sleep_score_sensitivity.csv", ["outcome"])
    compare_csv(required[4], ROOT / "expected_results/missingness_gee.csv", ["term"])
    compare_csv(required[5], ROOT / "expected_results/wellness_availability_ordinal_gee.csv", ["outcome"])


def main() -> None:
    validate_data()
    validate_privacy_scan()
    validate_outputs()
    print("Release validation passed: privacy guards, data invariants and key numerical outputs agree.")


if __name__ == "__main__":
    main()
