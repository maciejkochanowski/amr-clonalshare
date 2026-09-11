"""jsonio.py — strict, portable serialisation of the result dictionary.

``json.dumps`` will happily emit bare ``NaN`` and ``Infinity`` tokens, which
RFC 8259 forbids: JavaScript's ``JSON.parse``, R's ``jsonlite``, Go's
``encoding/json`` and Rust's ``serde_json`` all reject such a file. A results
artifact that only Python can read is not an interchange format.

``default=str`` is the other trap: it silently turns any numpy scalar or pandas
object that reaches the encoder into a quoted string, so a number becomes
``"0.42"`` with no error anywhere.

:func:`dumps` converts numpy and pandas types properly, maps non-finite floats
to ``null``, and passes ``allow_nan=False`` so that anything it missed fails
loudly instead of producing an unparseable file.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

__all__ = ["to_jsonable", "dumps", "write_json"]


def to_jsonable(obj: Any) -> Any:
    """Recursively convert numpy/pandas types and non-finite floats."""
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, np.ndarray):
        return [to_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (pd.Timestamp,)):
        return obj.isoformat()
    if isinstance(obj, pd.Series):
        return to_jsonable(obj.to_dict())
    if isinstance(obj, pd.DataFrame):
        return to_jsonable(obj.to_dict("records"))
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(
        f"cannot serialise {type(obj).__name__} to JSON; add a rule in "
        f"amr_clonalshare.jsonio.to_jsonable rather than stringifying it")


def dumps(obj: Any, *, indent: int = 2) -> str:
    """Serialise strictly: no NaN/Infinity tokens, no silent stringification."""
    return json.dumps(to_jsonable(obj), indent=indent, allow_nan=False,
                      ensure_ascii=False)


def write_json(obj: Any, path, *, indent: int = 2) -> Path:
    p = Path(path)
    p.write_text(dumps(obj, indent=indent), encoding="utf-8")
    return p
