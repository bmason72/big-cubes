"""Shared run-directory and config-metadata helpers for CLI workflows."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any


def _alpha_suffix(index: int) -> str:
    """0 -> a, 25 -> z, 26 -> aa, ..."""
    if index < 0:
        raise ValueError("index must be non-negative")
    chars = []
    n = index
    while True:
        n, rem = divmod(n, 26)
        chars.append(chr(ord("a") + rem))
        if n == 0:
            break
        n -= 1
    return "".join(reversed(chars))


def _date_tag(now: datetime | None = None) -> str:
    stamp = now or datetime.now()
    return f"{stamp.day}{stamp.strftime('%b').lower()}{stamp.strftime('%y')}"


def next_run_dir(
    pipeline_name: str,
    *,
    out_root: str | Path = "out",
    validation: bool = False,
    now: datetime | None = None,
) -> Path:
    """Return the next available run directory path."""
    root = Path(out_root)
    prefix = f"{'valid_' if validation else ''}{pipeline_name}_{_date_tag(now)}"
    for idx in range(26 * 27):
        candidate = root / f"{prefix}{_alpha_suffix(idx)}"
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"could not find free run dir for prefix {prefix}")


def _json_safe(obj: Any) -> Any:
    if is_dataclass(obj):
        return _json_safe(asdict(obj))
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, float) and not math.isfinite(obj):
        if math.isinf(obj):
            return "Infinity" if obj > 0 else "-Infinity"
        return "NaN"
    if callable(obj):
        return getattr(obj, "__name__", repr(obj))
    try:
        json.dumps(obj, allow_nan=False)
        return obj
    except TypeError:
        return repr(obj)
    except ValueError:
        return repr(obj)


def write_config_json(path: str | Path, payload: dict[str, Any]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(_json_safe(payload), indent=2) + "\n", encoding="utf-8")
