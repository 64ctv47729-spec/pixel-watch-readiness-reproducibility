#!/usr/bin/env bash
set -euo pipefail

package_root="$(cd "$(dirname "$0")" && pwd)"
python_bin="${HRV_PYTHON:-python3}"
export PYTHONPATH="$package_root/code${PYTHONPATH:+:$PYTHONPATH}"

"$python_bin" "$package_root/code/validate_release.py"
if [[ -f "$package_root/MANIFEST.sha256" ]]; then
  (cd "$package_root" && shasum -a 256 -c MANIFEST.sha256)
fi

