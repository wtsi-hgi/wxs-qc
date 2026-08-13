#!/usr/bin/env python3

"""
Shared helper to parse scripts/testdata-create/testdata.consts.

testdata.consts is a plain bash-sourceable file of `KEY="value"` lines, some of
which use `$(cd ... && pwd)` command substitution to derive absolute paths
relative to the file's own location (e.g. WORKDIR). To resolve those paths
correctly we actually `source` the file in a real bash subprocess rather than
trying to reimplement shell command substitution in Python.

No CLI here - just functions for the other scripts in this folder to import.
Standard library only (shlex, subprocess, pathlib).
"""

import shlex
import subprocess
from pathlib import Path

CONSTS_PATH = Path(__file__).with_name("testdata.consts")

# Keys that hold single-line, colon/space-delimited lists rather than plain
# scalar values; excluded from get_consts() and exposed via get_samples()/get_trios().
_LIST_KEYS = {"SAMPLES", "TRIOS"}


def _discover_keys(consts_path: Path) -> list[str]:
    """Find the KEY names defined in testdata.consts, in file order.

    Uses shlex only to keep the discovery structured (avoid ad hoc regex on the
    right-hand side); the actual values are resolved separately by sourcing the
    file in bash, since some of them contain `$(...)` command substitution that
    plain string parsing cannot evaluate.
    """
    keys = []
    for line in consts_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, rest = line.partition("=")
        key = key.strip()
        if not key.isidentifier() or not key.isupper():
            continue
        # Validate the right-hand side is a well-formed shell-quoted token
        # (raises if malformed), without trying to interpret its contents.
        shlex.split(rest)
        keys.append(key)
    return keys


def _resolve_values(consts_path: Path, keys: list[str]) -> dict[str, str]:
    """Source testdata.consts in bash and capture the resolved value of each key."""
    record_sep = "\x1e"
    script_lines = [f'source "{consts_path}"']
    script_lines += [f'printf "%s{record_sep}" "${{{key}}}"' for key in keys]
    result = subprocess.run(
        ["bash", "-c", "\n".join(script_lines)],
        check=True,
        capture_output=True,
        text=True,
    )
    values = result.stdout.split(record_sep)
    return dict(zip(keys, values, strict=False))


def _load_all() -> dict[str, str]:
    keys = _discover_keys(CONSTS_PATH)
    return _resolve_values(CONSTS_PATH, keys)


def get_consts() -> dict[str, str]:
    """Return all scalar (non-list) constants from testdata.consts as a dict."""
    return {key: value for key, value in _load_all().items() if key not in _LIST_KEYS}


def get_samples() -> list[dict[str, str]]:
    """Return the fixed sample list, parsed from the SAMPLES constant.

    Each entry is a dict with keys `sample_id`, `population`, `sex`, in the
    order they appear in testdata.consts.
    """
    raw = _load_all()["SAMPLES"]
    samples = []
    for token in raw.split():
        sample_id, population, sex = token.split(":")
        samples.append({"sample_id": sample_id, "population": population, "sex": sex})
    return samples


def get_trios() -> list[dict[str, str]]:
    """Return the fixed trio list, parsed from the TRIOS constant.

    Each entry is a dict with keys `family_id`, `child`, `father`, `mother`, in
    the order they appear in testdata.consts.
    """
    raw = _load_all()["TRIOS"]
    trios = []
    for token in raw.split():
        family_id, child, father, mother = token.split(":")
        trios.append({"family_id": family_id, "child": child, "father": father, "mother": mother})
    return trios
