#!/usr/bin/env python3
"""Write a SHA-256 manifest for the release package."""
from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def main() -> None:
    files = []
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if not path.is_file() or path.is_symlink():
            continue
        if relative.as_posix() == "MANIFEST.sha256":
            continue
        if any(part in {".venv", "__pycache__"} for part in relative.parts):
            continue
        files.append(path)
    lines = [f"{digest(path)}  {path.relative_to(ROOT).as_posix()}" for path in sorted(files)]
    (ROOT / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
    print(f"Manifested {len(lines)} files")


if __name__ == "__main__":
    main()

