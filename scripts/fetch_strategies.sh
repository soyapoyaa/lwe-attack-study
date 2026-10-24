#!/usr/bin/env bash
# fplll's BKZ pruning/preprocessing strategies (~8.6 MB).
# Needed only if `python3 src/strategies.py` reports "none".
# Source: https://github.com/fplll/fplll (LGPL-2.1+), file strategies/default.json,
# pinned to fplll commit b477504, the last change to that file (2020-09-30).
# The recorded runs used exactly this file; the checksum below verifies it.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data
curl -fSL -o data/bkz_strategies.json \
  https://raw.githubusercontent.com/fplll/fplll/b477504c88caa2ffb446f5b9b4b8af74695f468a/strategies/default.json
echo "f516b0a6f0c580cff72e1e2c3562c44dc6f17e8f99613e9e4020e35481b27a18  data/bkz_strategies.json" \
  | sha256sum -c -
python3 src/strategies.py
