#!/usr/bin/env bash
# fplll's BKZ pruning/preprocessing strategies (~8.6 MB).
# Needed only if `python3 src/strategies.py` reports "none".
# Source: https://github.com/fplll/fplll (LGPL-2.1+), file strategies/default.json
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data
curl -fSL -o data/bkz_strategies.json \
  https://raw.githubusercontent.com/fplll/fplll/master/strategies/default.json
python3 src/strategies.py
