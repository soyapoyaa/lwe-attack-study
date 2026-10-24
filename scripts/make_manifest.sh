#!/usr/bin/env bash
# Rebuild MANIFEST.sha256: checksums of the code, data and analysis outputs.
# Run from the repository root after any change to them, and commit the result.
# Check with: sha256sum -c MANIFEST.sha256
set -euo pipefail
cd "$(dirname "$0")/.."
git ls-files src analysis results scripts requirements.txt \
  | grep -v '^scripts/make_manifest.sh$' \
  | LC_ALL=C sort \
  | xargs sha256sum > MANIFEST.sha256
echo "MANIFEST.sha256: $(wc -l < MANIFEST.sha256) files"
