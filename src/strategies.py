"""
Locate fplll's BKZ strategies file (pruning + preprocessing parameters).

pip-installed fpylll wheels ship with `BKZ.DEFAULT_STRATEGY` pointing at the
*build machine's* path, which does not exist on your machine, so
`BKZ.Param(strategies=BKZ.DEFAULT_STRATEGY)` raises

    RuntimeError: Cannot open strategies file.

Resolution order:
  1. $BKZ_STRATEGIES
  2. BKZ.DEFAULT_STRATEGY, if that file actually exists (conda / Sage / source builds)
  3. <repo>/data/bkz_strategies.json   (run scripts/fetch_strategies.sh once)
  4. None -- BKZ still runs, but without preprocessing or pruning, so it is
     several times slower. A warning is printed once.

Running without strategies is a *performance* difference, not a correctness
one, but it changes measured wall-clock times, so never mix strategies and
no-strategies runs in one table. `describe()` records which was used.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

from fpylll import BKZ, load_strategies_json

_REPO_DATA = Path(__file__).resolve().parent.parent / "data" / "bkz_strategies.json"
_cache: dict = {}


def strategies_path() -> Path | None:
    env = os.environ.get("BKZ_STRATEGIES")
    if env and Path(env).exists():
        return Path(env)
    try:
        default = BKZ.DEFAULT_STRATEGY
        default = default.decode() if isinstance(default, bytes) else default
        if default and Path(default).exists():
            return Path(default)
    except Exception:
        pass
    if _REPO_DATA.exists():
        return _REPO_DATA
    return None


def get_strategies():
    """Loaded strategy vector, or None."""
    if "value" in _cache:
        return _cache["value"]
    path = strategies_path()
    if path is None:
        warnings.warn(
            "No BKZ strategies file found: running BKZ without preprocessing or "
            "pruning (correct but slow). Run scripts/fetch_strategies.sh once.",
            RuntimeWarning, stacklevel=2,
        )
        _cache["value"] = None
        _cache["path"] = None
    else:
        _cache["value"] = load_strategies_json(str(path).encode())
        _cache["path"] = str(path)
    return _cache["value"]


def describe() -> str:
    get_strategies()
    return _cache.get("path") or "none (no preprocessing/pruning)"


def bkz_param(block_size: int, **kwds) -> BKZ.Param:
    st = get_strategies()
    flags = kwds.pop("flags", BKZ.GH_BND | BKZ.AUTO_ABORT)
    if st is None:
        return BKZ.Param(block_size=block_size, flags=flags, **kwds)
    return BKZ.Param(block_size=block_size, strategies=st, flags=flags, **kwds)


if __name__ == "__main__":
    print("strategies:", describe())
    print("param:", bkz_param(30))
