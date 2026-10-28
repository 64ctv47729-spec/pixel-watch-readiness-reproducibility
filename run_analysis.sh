#!/usr/bin/env bash
set -euo pipefail

package_root="$(cd "$(dirname "$0")" && pwd)"
python_bin="${HRV_PYTHON:-python3}"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/hrv_readiness_matplotlib"
export PYTHONPATH="$package_root/code${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "$package_root/output/audit" "$package_root/output/models" \
  "$package_root/output/reviewer" "$package_root/results/corrected/aim3"
cp "$package_root/data/analytic_dataset.csv" "$package_root/output/analytic_dataset.csv"
cp "$package_root/data/duplicate_questionnaire_dates.csv" \
  "$package_root/output/audit/duplicate_questionnaire_adjudication.csv"
cp "$package_root/data/audit_summary.json" "$package_root/output/audit/audit_summary.json"

"$python_bin" "$package_root/code/independent_models.py"
"$python_bin" "$package_root/code/missingness.py"
"$python_bin" "$package_root/code/missingness_ipw_sensitivity.py"
"$python_bin" "$package_root/code/timing_window_sensitivity.py"
"$python_bin" "$package_root/code/reviewer_additional_analyses.py"
"$python_bin" "$package_root/code/aim3_corrected.py"
"$python_bin" "$package_root/code/aim3_sleep_score_sensitivity.py"
"$python_bin" "$package_root/code/make_corrected_tables_figures.py"
"$python_bin" "$package_root/code/validate_release.py"
"$python_bin" "$package_root/code/make_manifest.py"

echo "Reproduction and validation completed: $package_root"

