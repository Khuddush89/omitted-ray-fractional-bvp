#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
command -v pdftoppm >/dev/null || { echo "Install Poppler first: brew install poppler"; exit 1; }
mkdir -p figures/previews
for figure in figures/fig*.pdf; do
  stem="$(basename "$figure" .pdf)"
  pdftoppm -png -singlefile -r 180 "$figure" "figures/previews/$stem"
done
echo "All figure previews updated."
